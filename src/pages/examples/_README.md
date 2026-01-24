# Nextflow Examples

This directory contains example Nextflow pipelines with their documentation pages.

## Structure

Each example follows a consistent structure:

```
example-name/
├── index.mdx           # Documentation page with explanation
├── _main.nf            # Nextflow pipeline script
├── _nextflow_run_output.log  # Captured pipeline execution output
└── data/               # Input data files (if needed)
    └── sample.fa       # Sample input file
```

## File Naming Convention

- **`_main.nf`**: Pipeline script (prefixed with `_` to avoid conflicts with actual `main.nf`)
- **`_nextflow_run_output.log`**: Raw terminal output from running the pipeline
- **`index.mdx`**: Documentation page that imports and displays both the pipeline code and execution output

## Adding New Examples

1. Create a new directory with a descriptive name
2. Add the pipeline script as `_main.nf`
3. Create any necessary input data in a `data/` subdirectory
4. Run the pipeline and capture output with ANSI colors:
   ```bash
   # For colorized output that renders properly with lang="ansi"
   nextflow run _main.nf > _nextflow_run_output.log 2>&1

   # Then convert text escape sequences to actual binary escape characters
   # Replace \x1B[ patterns with actual escape characters using printf
   printf "nextflow run _main.nf\n\n\x1B[1;42m N E X T F L O W \x1B[0m  ~  version X.X.X\n..." > _nextflow_run_output.log
   ```
5. Create an `index.mdx` file that imports both files:
   ```javascript
   import pipelineCode from "./_main.nf?raw";
   import terminalOutput from "./_nextflow_run_output.log?raw";
   ```

## Documentation Pages

Each `index.mdx` file should include:

- Clear explanation of what the pipeline does
- Key concepts demonstrated
- Code blocks showing both the pipeline and execution output
- Usage instructions

The pages use Expressive Code for syntax highlighting and the `ExampleLayout` for consistent styling.

## ANSI Color Support

Terminal output uses `lang="ansi"` to render colorized output. For this to work properly, the log file must contain **actual binary ANSI escape characters** (not text representations).

### Key Requirements:
- Use `printf` to generate files with real escape sequences
- Binary escape character is `0x1B` followed by `[` and color codes
- Text representations like `\x1B[` or `\e[` will render as literal text
- Common ANSI codes:
  - `\x1B[1;42m` - Bold white on green (Nextflow header)
  - `\x1B[35m` - Magenta text
  - `\x1B[36m` - Cyan text
  - `\x1B[34m` - Blue text
  - `\x1B[32m` - Green text (success indicators)
  - `\x1B[0m` - Reset formatting