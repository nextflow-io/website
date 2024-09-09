---
title: Machine Learning pipeline
layout: "@layouts/MarkdownPage.astro"
---

<div class="blg-summary example">
<h3>Machine Learning pipeline</h3>

<p class="text-muted">
    This example shows how to put together a basic Machine Learning pipeline.
</p>

```groovy
#!/usr/bin/env nextflow

params.dataset_name = 'wdbc'
params.train_models = ['dummy', 'gb', 'lr', 'mlp', 'rf']
params.outdir = 'results'

workflow {
    // fetch dataset from OpenML
    ch_datasets = fetch_dataset(params.dataset_name)

    // split dataset into train/test sets
    (ch_train_datasets, ch_predict_datasets) = split_train_test(ch_datasets)

    // perform training
    (ch_models, ch_train_logs) = train(ch_train_datasets, params.train_models)

    // perform inference
    ch_predict_inputs = ch_models.combine(ch_predict_datasets, by: 0)
    (ch_scores, ch_predict_logs) = predict(ch_predict_inputs)

    // select the best model based on inference score
    ch_scores
        | max {
            new JsonSlurper().parse(it[2])['value']
        }
        | subscribe { dataset_name, model_type, score_file ->
            def score = new JsonSlurper().parse(score_file)
            println "The best model for ${dataset_name} was ${model_type}, with ${score['name']} = ${score['value']}"
        }
}

// view the entire code on GitHub ...

```

</div>

### Synopsis

This example shows how to put together a basic Machine Learning pipeline. It fetches a dataset from OpenML, trains a variety of machine learning models on a prediction target, and selects the best model based on some evaluation criteria.

### Try it

This pipeline is available on the [nextflow-io/ml-hyperopt](https://github.com/nextflow-io/ml-hyperopt) GitHub repository.

An active internet connection and Docker are required for Nextflow to download the pipeline and the necessary Docker images to run the pipeline within containers. The data used by this pipeline is stored on the GitHub repository and will download automatically.

To try this pipeline:

1. Follow the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html#install-nextflow) to install Nextflow.
2. Follow the [Docker installation guide](https://docs.docker.com/get-started/get-docker/) to install Docker.
3. Launch the pipeline:

   nextflow run nextflow-io/ml-hyperopt -profile wave

**NOTE**: Nextflow 22.10.0 or newer is required to run this pipeline with Wave.
