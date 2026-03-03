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