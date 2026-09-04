# Negative phonon branch study

This folder is independent of `TDEP_300K/codes-workflow`. It snapshots the
currently prepared TDEP trajectory into multiple cutoff, stride, time-window,
and cumulative-sampling fits, then creates comparison plots and `REPORT.md`.

Run again with:

```bash
/opt/anaconda3/bin/python TDEP_negative_branch_study/run_study.py
```

The study intentionally tests cutoffs only through 5.5 Å because this
supercell's maximum valid pair cutoff is approximately 5.67 Å.
