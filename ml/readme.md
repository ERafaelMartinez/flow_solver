# SimML

Automatically train multiple neural networks against the LidCavity scenario and analyse them to choose the best one.

## Getting started

1. Make sure you are using `python v3.11`.
2. Set up your virtual environment and install all dependencies in `requirements.txt`.
3. Call `python main.py` to train all existing models and create an analysis under `assets/runs/run-{timestamp}`

## How to add a new network/model

To add a new network, simple create a new class in `networks/` that inherits from `BaseModel`. 
The network will then be automatically detected, trained and analysed the next time you call `python main.py`

> **tip** to start quickly you can simply duplicate `ReferenceCNN.py` 
> but make sure to change the name of both the file and the class

> **tip** you can exclude a model from being trained and tested by simply prefixing the name 
> of the file with an `x`. So `xReferenceCNN` will be ignored and won't be trained. To train this network
> simply remove the x


## How to change the hyperparameters

Not yet implemented
<!-- TODO: Still not possible -->

## How to change the training method

All training is being done by the `trainer.py` module.

## How to change the analysis

The analysis consists of 
1. training analysis

The training analysis is being done while training the model. After the `Trainer` is done training a model
it reports back the training loss and validation loss per epoch.
These metrics can be changed and modified in the `trainer.py` module. Make sure to always return a `TrainResult` object
with the metrics needed

2. testing analysis

The testing analysis is done after the model is finished training against the testing dataset. 
After the `Tester` is done testing the model, it currently only reports some metrics. 
These metrics can be changed and modified in the `tester.py` module. Make sure to always return a `TestResult` object
with the metrics needed


After the `Trainer` and the `Tester` are finished, the `Reporter` takes all the metrics and generates reports under the
run folder (e.g. `assets/runs/run-{timestamp}`)
