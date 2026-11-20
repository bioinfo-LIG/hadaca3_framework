# Contributing

Pull requests are welcome. Please report bugs and propose features in the
[issues](https://github.com/bioinfo-LIG/hadaca3_framework/issues).

## Adding a method

1. **Declare the block** in the yml file of its module (`function_metadata_and_selection/template_complete_run/`).

   ```yaml
   myMethod:
     path: function_blocks/<module>/my_method.R
     short_name: MyMth      # 4-5 characters, unique
     language: R            # or python (early integration only)
     omic: [RNA, MET]       # TO CHECK against the current README
   ```

   - **Required fields:** `path`, `short_name`, `omic`.
   - **Optional fields:** `dependency`, `create`, `need` / `omic_need`.

2. **Implement it** in `function_blocks/<module>/`, with the signature of its module. For example, for a
   pre-processing block: `program_block_PP(data, path_og_dataset='', omic='')`. `ppID` (identity) is a template.
3. **Document it.**
   - Add a line to [function_blocks/readme.md](function_blocks/readme.md).
   - If it changes or replaces an existing block, add a note at the top of the old block and a line to
     [KNOWN_ISSUES.md](KNOWN_ISSUES.md).
   - **Do not delete a block used in a published analysis**: deprecate it, so that the published results stay
     reproducible.

## Continuous integration

- **On every commit:** a partial run (`-resume`).
- **Daily at 3 am:** a full run, if the repository changed.
- **Incremental run:** it uses `function_metadata_and_selection/CI_functions_selection`.
- **Skipping the CI:** put `[skip ci]`, `[ci skip]`, `[no ci]`, `[skip actions]` or `[actions skip]` in the commit
  message.
