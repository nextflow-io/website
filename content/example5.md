title=Machine Learning pipeline
date=2022-10-31
type=page
status=published
syntaxhighlighter=yes
~~~~~~

<div class="blg-summary example">
<h3><a href="javascript:void(0)">Machine Learning pipeline</a></h3>

<p class="text-muted">
    This example shows how to put together a basic Machine Learning pipeline. It fetches a dataset from OpenML, trains a variety of machine learning models on a prediction target, and selects the best model based on some evaluation criteria.
</p>

<script type="syntaxhighlighter" class="brush: groovy">
<![CDATA[
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

]]>
</script>
</div>


### Try it in your computer

To run this pipeline on your computer, you will need:

* Unix-like operating system
* Java 11 (or higher)
* Docker

Install Nextflow by entering the following command in the terminal:

    $ curl -fsSL get.nextflow.io | bash

Then launch the pipeline with this command:

    $ nextflow run ml-hyperopt -profile wave

It will automatically download the pipeline [GitHub repository](https://github.com/nextflow-io/ml-hyperopt) and build a Docker image on-the-fly using [Wave](https://seqera.io/wave/), thus the first execution may take a few minutes to complete depending on your network connection.

__NOTE__: Nextflow 22.10.0 or newer is required to run this pipeline with Wave.
