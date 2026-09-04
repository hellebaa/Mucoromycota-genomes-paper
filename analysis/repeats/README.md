# Repeat Annotation Summaries

This folder contains scripts for summarizing repeat annotations used in the
Mucoromycota comparative genomics paper.

## Figure 1 De Novo Repeat Summary

The draft Figure 1 repeat panel used summary statistics parsed from
RepeatMasker `.tbl` files generated during the species-specific de novo repeat
annotation runs.

Historical Saga input directory:

```text
/cluster/projects/nn8013k/results/helle/repeats/repeat_tbl
```

Historical output:

```text
repeat_summary_table.tsv
```

Reusable command:

```bash
python analysis/repeats/scripts/parse_repeatmasker_tbl.py \
  --input-dir /cluster/projects/nn8013k/results/helle/repeats/repeat_tbl \
  --output repeat_summary_table.tsv
```

When `MucCirc3` de novo RepeatModeler/RepeatMasker output is complete and on
NeLS/Saga, rerun this parser so the Figure 1 repeat table includes `MucCirc3`.
