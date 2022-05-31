
chromosome = ""

with open("results/qc/included_regions/TAIR10.col0.bg") as bedgraph:
    for line in bedgraph:
        chrom, start, end, value = line.strip().split("\t")

        if chrom != chromosome:
            print(f"fixedStep chrom={chrom} start=1 step=1 span=1")
            chromosome = chrom
        
        for i in range(int(start), int(end)):
            print(value)
