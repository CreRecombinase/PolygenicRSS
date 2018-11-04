with open(snakemake.output[0],"w") as out:
    for f in snakemake.params.in_pref[1:]:
        out.write(str(f)+"\n")
