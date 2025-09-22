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
4. Run the pipeline and capture output:
   ```bash
   nextflow run _main.nf > _nextflow_run_output.log 2>&1
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