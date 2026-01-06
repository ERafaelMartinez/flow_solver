import logging
import torch
from hyperparameters import Hyperparameters
from normalizer import Normalizer

logger = logging.getLogger(__name__)


def load_data(normalizer: Normalizer) -> tuple[torch.utils.data.DataLoader,
                                               torch.utils.data.DataLoader,
                                               torch.utils.data.DataLoader]:
    n_samples = 101
    train_ratio = 0.8
    validation_ratio = 0.1
    batch_size = 10

    inputs = torch.load("data/inputs.pt")
    outputs = torch.load("data/outputs.pt")

    (inputs_normalized, outputs_normalized) = normalizer.normalize(inputs, outputs)

    # We need to make sure that we do not use more samples than what we already
    # have unless we want to augment the data somehow.
    n_samples = min(n_samples, inputs.shape[0])
    n_train = int(train_ratio * n_samples)
    n_validation = int(validation_ratio * n_samples)
    n_test = n_samples - n_train - n_validation

    logger.info(f"Splitting {n_samples} samples into: {n_train} training, {n_validation} validation, {n_test} testing")

    # Generate a random permutation of indices from 0 to n_samples-1
    indices = torch.randperm(inputs.shape[0])

    # Slice the indices for train and test sets
    train_indices = indices[:n_train]
    validation_indices = indices[n_train:(n_train+n_validation)]
    testing_indices = indices[(n_train+n_validation):]

    logger.debug(f"Indices used for training\n{train_indices}")
    logger.debug(f"Indices used for validation\n{validation_indices}")
    logger.debug(f"Indices used for testing\n{testing_indices}")

    # create tensor dataset and dataloader
    train_dataset = torch.utils.data.TensorDataset(
        inputs_normalized[train_indices], outputs_normalized[train_indices])
    training_dataloader = torch.utils.data.DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True)

    # We do not want to shuffle the validation dataset
    validation_dataset = torch.utils.data.TensorDataset(
        inputs_normalized[validation_indices], outputs_normalized[validation_indices])
    validation_dataloader = torch.utils.data.DataLoader(
        validation_dataset,
        shuffle=False)

    # We do not want to shuffle the testing dataset
    testing_dataset = torch.utils.data.TensorDataset(
        inputs_normalized[testing_indices], outputs_normalized[testing_indices])
    testing_dataloader = torch.utils.data.DataLoader(
        testing_dataset,
        shuffle=False)

    return (training_dataloader, validation_dataloader, testing_dataloader)
