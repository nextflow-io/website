#!/usr/bin/env nextflow

params.input_data = "${baseDir}/data/dataset.csv"
params.test_size = 0.3
params.models = ['random_forest', 'svm', 'logistic_regression']

workflow {
    // Create input channel from dataset
    ch_dataset = Channel.fromPath(params.input_data)

    // Split dataset into training and test sets
    split_dataset(ch_dataset)

    // Train multiple models in parallel
    ch_models_input = Channel.of(params.models)
    train_models(split_dataset.out.train, ch_models_input)

    // Evaluate each model on test data
    evaluate_models(train_models.out, split_dataset.out.test)

    // Find the best performing model
    evaluate_models.out.view { "Model: ${it[0]}, Accuracy: ${it[1]}" }
}

process split_dataset {
    input:
    path dataset

    output:
    path 'train_data.csv', emit: train
    path 'test_data.csv', emit: test

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from sklearn.model_selection import train_test_split

    # Load dataset
    data = pd.read_csv('${dataset}')
    X = data.iloc[:, :-1]  # Features
    y = data.iloc[:, -1]   # Target

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=${params.test_size}, random_state=42
    )

    # Save splits
    train_data = pd.concat([X_train, y_train], axis=1)
    test_data = pd.concat([X_test, y_test], axis=1)

    train_data.to_csv('train_data.csv', index=False)
    test_data.to_csv('test_data.csv', index=False)
    """
}

process train_models {
    input:
    path train_data
    each model_type

    output:
    tuple val(model_type), path("${model_type}_model.pkl")

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import pickle
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.linear_model import LogisticRegression

    # Load training data
    data = pd.read_csv('${train_data}')
    X = data.iloc[:, :-1]
    y = data.iloc[:, -1]

    # Select and train model
    if '${model_type}' == 'random_forest':
        model = RandomForestClassifier(random_state=42)
    elif '${model_type}' == 'svm':
        model = SVC(random_state=42)
    elif '${model_type}' == 'logistic_regression':
        model = LogisticRegression(random_state=42)

    # Train the model
    model.fit(X, y)

    # Save the model
    with open('${model_type}_model.pkl', 'wb') as f:
        pickle.dump(model, f)
    """
}

process evaluate_models {
    input:
    tuple val(model_type), path(model_file)
    path test_data

    output:
    tuple val(model_type), stdout

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import pickle
    from sklearn.metrics import accuracy_score

    # Load test data
    data = pd.read_csv('${test_data}')
    X_test = data.iloc[:, :-1]
    y_test = data.iloc[:, -1]

    # Load and evaluate model
    with open('${model_file}', 'rb') as f:
        model = pickle.load(f)

    # Make predictions and calculate accuracy
    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)

    print(f"{accuracy:.4f}")
    """
}
