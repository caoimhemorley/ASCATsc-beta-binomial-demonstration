#!/usr/bin/env python3

import argparse

import csv

import gzip

import os

import re

import subprocess

from pathlib import Path

def opener(path):

    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def locate_one(root, patterns, label):

    matches = sorted({

        p for pattern in patterns

        for p in root.glob(pattern)

        if p.is_file()

    })

    if len(matches) != 1:

        raise RuntimeError(

            f"{label}: expected exactly one file; found {len(matches)}: {matches[:8]}"

        )

    return matches[0]

def locate_ref(chrom, container, reference_dir):

    cmd = [

        "apptainer", "exec", container,

        "sh", "-c", f"find {reference_dir} -maxdepth 2 -type f"

    ]

    output = subprocess.check_output(cmd, text=True)

    files = [s.strip() for s in output.splitlines()]

    matches = [

        f for f in files

        if re.match(rf"^(?:chr)?{chrom}(?:[._-]|$)", Path(f).name, re.I)

        and f.endswith((".vcf.gz", ".bref3", ".vcf"))

    ]

    if len(matches) != 1:

        raise RuntimeError(

            f"Reference for chr{chrom} ambiguous/missing: {matches}"

        )

    return matches[0]

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--donor", required=True, help="Sample identifier")

    parser.add_argument("--sex", choices=["male", "female"], default="female")
    parser.add_argument("--chrom", required=True, choices=[str(x) for x in range(1, 23)] + ["X"])

    parser.add_argument("--min-depth", type=int, default=8)

    parser.add_argument("--min-alt-reads", type=int, default=3)

    parser.add_argument("--min-het-fraction", type=float, default=0.20)

    parser.add_argument("--max-hom-minor-reads", type=int, default=1)

    parser.add_argument("--max-hom-minor-fraction", type=float, default=0.05)

    parser.add_argument("--threads", type=int, default=4)

    parser.add_argument(

        "--out-root",

        default=None

    )

    parser.add_argument("--input-root", required=True, help="Sample-specific directory containing bulk allele counts")
    parser.add_argument("--loci-dir", required=True, help="Directory of per-chromosome GRCh38 loci tables")
    parser.add_argument("--beagle-container", required=True, help="Apptainer SIF with Beagle, reference panel and maps")
    parser.add_argument("--reference-dir", default="/home/ref_files/hg38")
    parser.add_argument("--map-dir", default="/home/mapfiles")
    parser.add_argument("--beagle-jar", default="/home/beagle.22Jul22.46e.jar")
    parser.add_argument("--prepare-only", action="store_true")

    args = parser.parse_args()

    if not 0 < args.min_het_fraction < 0.5:

        raise ValueError("Invalid heterozygote fraction")

    root = Path(args.input_root)


    allele_dir = root / "allele_count"

    loci_candidates = [
        file for file in Path(args.loci_dir).iterdir()
        if file.is_file()
        and re.search(rf"(^|[._-])chr{args.chrom}([._-]|$)", file.name)
        and (file.name.endswith(".txt") or file.name.endswith(".txt.gz"))
    ]
    if len(loci_candidates) != 1:
        raise RuntimeError(
            f"GRCh38 loci table chr{args.chrom}: expected exactly one; found "
            f"{len(loci_candidates)}: {sorted(loci_candidates)}"
        )
    loci_file = loci_candidates[0]

    allele_file = locate_one(

        allele_dir,

        [

            f"*_chr{args.chrom}.txt",

            f"*_chr{args.chrom}.txt.gz"

        ],

        "Bulk allele counts"

    )

    out_dir = (Path(args.out_root) if args.out_root else root / "beagle_confident_v1") / f"chr{args.chrom}"

    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = out_dir / f"{args.donor}_chr{args.chrom}"

    target = Path(str(prefix) + ".confident_input.vcf")

    summary = Path(str(prefix) + ".summary.tsv")

    final = Path(str(prefix) + ".beagleoutput.vcf.gz")

    final_index = Path(str(final) + ".tbi")

    if final.is_file() and final.stat().st_size > 0 and final_index.is_file() and final_index.stat().st_size > 0:

        check = subprocess.run(

            ["tabix", "-l", str(final)],

            capture_output=True,

            text=True

        )

        if check.returncode == 0 and check.stdout.strip():

            print(

                f"SKIPPING completed {args.donor} chr{args.chrom}: {final}",

                flush=True

            )

            return

        print(

            f"Existing output invalid; rerunning chr{args.chrom}",

            flush=True

        )

    counts = {}

    with opener(allele_file) as handle:

        reader = csv.DictReader(handle, delimiter="\t")

        required = {"POS", "Count_A", "Count_C", "Count_G", "Count_T"}

        if not required.issubset(reader.fieldnames or []):

            raise ValueError(

                f"Unexpected allele count columns: {reader.fieldnames}"

            )

        for row in reader:

            pos = int(row["POS"])

            if pos in counts:

                raise ValueError(f"Duplicate allele count position: {pos}")

            counts[pos] = tuple(

                int(row[f"Count_{base}"]) for base in "ACGT"

            )

    stats = {

        "total_input": 0,

        "het": 0,

        "hom_ref": 0,

        "hom_alt": 0,

        "ambiguous": 0,

        "no_count": 0,

        "non_snv": 0,

        "depth_below_min": 0,

        "other_base_fraction_gt_05": 0

    }

    accepted_hets = set()

    accepted_homs = {}

    is_male_x = args.chrom == "X" and args.sex == "male"

    with opener(loci_file) as source, open(target, "w") as dest:
        dest.write("##fileformat=VCFv4.2\n")
        dest.write("##reference=GRCh38\n")
        dest.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        dest.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.donor + "\n")
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                raise ValueError(f"Expected chromosome, position, REF, ALT in {loci_file}")
            chrom, pos_text, ref, alt = fields[:4]
            if not pos_text.isdigit():
                if pos_text.lower() in ("pos", "position"):
                    continue
                raise ValueError(f"Invalid position in {loci_file}: {pos_text}")
            pos = int(pos_text)
            chrom = chrom.removeprefix("chr")
            if chrom != args.chrom:
                continue
            stats["total_input"] += 1
            if args.chrom == "X" and not (2781480 <= pos <= 155701382):
                continue
            ref, alt = ref.upper(), alt.upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                stats["non_snv"] += 1
                continue
            if pos not in counts:
                stats["no_count"] += 1
                continue
            base_counts = dict(zip("ACGT", counts[pos]))
            ref_count, alt_count = base_counts[ref], base_counts[alt]
            depth = ref_count + alt_count
            if depth < args.min_depth:
                stats["depth_below_min"] += 1
                continue
            other = sum(counts[pos]) - depth
            if other / (depth + other) > 0.05:
                stats["other_base_fraction_gt_05"] += 1
                continue
            minor = min(ref_count, alt_count)
            if minor >= args.min_alt_reads and minor / depth >= args.min_het_fraction:
                gt = "0/1"
                stats["het"] += 1
                accepted_hets.add((chrom, pos, ref, alt))
            elif minor <= args.max_hom_minor_reads and minor / depth <= args.max_hom_minor_fraction:
                gt = "0/0" if ref_count > alt_count else "1/1"
                stats["hom_ref" if ref_count > alt_count else "hom_alt"] += 1
                if is_male_x:
                    accepted_homs[(chrom, pos, ref, alt)] = gt
            else:
                stats["ambiguous"] += 1
                continue
            dest.write("\t".join((chrom, str(pos), ".", ref, alt, ".", "PASS", ".", "GT", gt)) + "\n")

    with open(summary, "w") as handle:

        handle.write("metric\tcount\n")

        for metric, value in stats.items():

            handle.write(f"{metric}\t{value}\n")

    print(f"Donor: {args.donor}", flush=True)

    print(f"Chromosome: {args.chrom}", flush=True)

    print(f"Allele counts: {allele_file}", flush=True)

    print(f"Loci table: {loci_file}", flush=True)

    print(f"Genotypes: {stats}", flush=True)

    if not is_male_x and stats["het"] < 1000:

        raise RuntimeError(

            "Fewer than 1000 confident heterozygotes; investigate data"

        )

    if args.prepare_only:

        print(f"Prepared VCF: {target}", flush=True)

        print(f"QC summary: {summary}", flush=True)

        return

    ref = locate_ref(args.chrom, args.beagle_container, args.reference_dir)

    genetic_map = f"{args.map_dir}/plink.chr{args.chrom}.GRCh38.map"

    output_prefix = str(prefix) + ".beagle"

    cmd = [

        "apptainer", "exec",

        args.beagle_container,

        "java", "-Xmx12g",

        "-jar", args.beagle_jar,

        f"gt={target}",

        f"ref={ref}",

        f"map={genetic_map}",

        f"out={output_prefix}",

        f"nthreads={args.threads}",

        "impute=false"

    ]

    print("Running:", " ".join(cmd), flush=True)

    subprocess.run(cmd, check=True)

    phased = Path(output_prefix + ".vcf.gz")

    if not phased.exists():

        raise FileNotFoundError(phased)

    final = Path(str(prefix) + ".beagleoutput.vcf.gz")

    final_plain = Path(str(prefix) + ".beagleoutput.vcf")

    retained = 0

    seen = 0

    with gzip.open(phased, "rt") as source, open(final_plain, "w") as dest:

        for line in source:

            if line.startswith("#"):

                dest.write(line)

                continue

            fields = line.rstrip("\n").split("\t")

            seen += 1

            key = (

                fields[0],

                int(fields[1]),

                fields[3],

                fields[4]

            )

            if is_male_x:

                if key not in accepted_homs:

                    continue

                fields[9] = accepted_homs[key].replace("/", "|")

            else:

                if key not in accepted_hets:

                    continue

                gt = fields[9].split(":", 1)[0]

                if gt not in ("0|1", "1|0"):

                    continue

            if fields[7] not in (".", ""):

                fields[7] = (

                    ";".join(

                        x for x in fields[7].split(";") if x != "IMP"

                    ) or "."

                )

            dest.write("\t".join(fields) + "\n")

            retained += 1

    print(

        f"Beagle records: {seen}; "

        f"confident phased heterozygotes retained: {retained}",

        flush=True

    )

    if not is_male_x and retained < 1000:

        raise RuntimeError(

            "Too few phased heterozygotes; investigate reference matching"

        )

    subprocess.run(

        ["bgzip", "-f", str(final_plain)],

        check=True

    )

    subprocess.run(

        ["tabix", "-f", "-p", "vcf", str(final)],

        check=True

    )

    print(f"Final phased VCF: {final}", flush=True)

    print(f"Completed {args.donor} chr{args.chrom}: {retained} retained records", flush=True)

if __name__ == "__main__":

    main()
