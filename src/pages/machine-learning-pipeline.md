---
title: Machine Learning pipeline
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>Machine Learning pipeline</h2>

<p class="">
    This pipeline shows how to put together a basic machine learning pipeline. It fetches datasets from OpenML, trains a variety of machine learning models on each dataset, evaluates each model against a test set, and reports the best model for each dataset. The excerpt below shows the core of the pipeline.
</p>

```nextflow
nextflow.enable.moduleBinaries = true
nextflow.enable.types = true

include { FETCH_DATASET } from './modules/fetch_dataset'
include { SPLIT_TRAIN_TEST } from './modules/split_train_test'
include { VISUALIZE as VISUALIZE_TRAIN } from './modules/visualize'
include { VISUALIZE as VISUALIZE_TEST } from './modules/visualize'
include { TRAIN } from './modules/train'
include { EVALUATE } from './modules/evaluate'

params {
    fetch_datasets: String?
    train_test_splits: Path?
    train_models: String?
    pretrained_models: Path?
}

workflow {
    main:
    // fetch and split datasets if specified
    if( params.fetch_datasets != null ) {
        ch_dataset_names = channel.fromList(params.fetch_datasets.tokenize(','))
        ch_datasets = FETCH_DATASET(ch_dataset_names)
        ch_train_test_splits = SPLIT_TRAIN_TEST(ch_datasets)
    }

    // ...

    // separate training and test data
    ch_train_datasets = ch_train_test_splits.map { r ->
        record(dataset_name: r.dataset_name, meta: r.meta, data: r.data_train)
    }
    ch_test_datasets = ch_train_test_splits.map { r ->
        record(dataset_name: r.dataset_name, meta: r.meta, data: r.data_test)
    }

    // ...

    // train new models if specified
    def ch_models: Channel<Model>
    if( params.train_models != null ) {
        model_types = params.train_models.tokenize(',')
        ch_train_inputs = ch_train_datasets.flatMap { r ->
            model_types.collect { model_type ->
                r + record(model_type: model_type)
            }
        }
        ch_models = TRAIN(ch_train_inputs)
    }

    // ...

    // evaluate each model against test dataset
    ch_evaluate_inputs = ch_models.join(ch_test_datasets, by: 'dataset_name')
    ch_evals = EVALUATE(ch_evaluate_inputs).map { r ->
        record(
            model_type: r.model_type,
            dataset_name: r.dataset_name,
            score_name: r.score.name,
            score: r.score.value,
            logs: r.logs,
        )
    }

    // report the best model for each dataset based on evaluation score
    ch_evals
        .map { r -> tuple(r.dataset_name, r) }
        .groupBy()
        .subscribe { dataset_name, evals ->
            def best = evals.max { r -> r.score }
            printf "The best model for dataset '${dataset_name}' was '${best.model_type}' (${best.score_name} = %.3f)\n", best.score
        }

    // ...
}

// ...
```

</div>

### Synopsis

The pipeline imports a process for each step from the `modules` directory:

- `FETCH_DATASET`: downloads a dataset from [OpenML](https://www.openml.org/)
- `SPLIT_TRAIN_TEST`: splits each dataset into a training set and a test set
- `VISUALIZE`: plots each training and test set
- `TRAIN`: trains a model of a given type (e.g. random forest or logistic regression) on a training set
- `EVALUATE`: scores each trained model against the corresponding test set

The `workflow` block uses records to keep each dataset name together with its files. It uses the `flatMap` operator to create a training task for every combination of dataset and model type, so that all models are trained in parallel, and the `join` operator to match each model with the test set for the same dataset. Finally, it groups the evaluation scores by dataset and reports the best model for each one.

You can also provide your own train/test splits or pre-trained models instead of fetching datasets and training new models.

### Get started

To run this pipeline:

1. [Install Nextflow](https://docs.seqera.io/nextflow/install) (version 26.10 or later).
2. [Install Docker](https://docs.docker.com/get-started/get-docker/).
3. Run the pipeline directly from its GitHub repository with the test data:

   ```bash
   nextflow run nextflow-io/hyperopt -profile test,wave
   ```

The `wave` profile uses [Wave](https://docs.seqera.io/nextflow/wave) to build a container from the pipeline's Conda environment on the fly, so the first run may take a few minutes. When the run completes, the datasets, trained models, and evaluation results are saved in the `results` directory.

See the [hyperopt](https://github.com/nextflow-io/hyperopt) GitHub repository for all of the pipeline code.
