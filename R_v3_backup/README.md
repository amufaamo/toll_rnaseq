# root and legacy archive

This directory contains retired v3 modules, root `app.R` snapshots, old R
session files, and superseded project-note copies. None of these files are
part of the v4 release.

Current entry points remain at `../app.R` for the integrated v3 app and
`../R/` for the v4 package. Archived files are retained only for manual
recovery and are ignored by Git, except for this README.

Archive layout:

- `module_*.R`: retired v3 modules collected on 2026-05-29.
- `app.R.bak_*`: dated or named snapshots formerly stored in the repository root.
- `root_session_*.RData`: old root-level R workspaces and temporary workspaces.
- `root_*.Rhistory`: old root-level command history.
- `root_Rplots_*.pdf`: unclassified plot output formerly stored in the root.
- `before_count_*.RData` and `before_count_*.Rhistory`: old session state
  removed from the two legacy preprocessor implementations.
- `repo_*.md`: superseded in-repository project notes. The canonical notes are
  `/mnt/g/マイドライブ/toll_project.md` and
  `/mnt/g/マイドライブ/toll_log.md`.
