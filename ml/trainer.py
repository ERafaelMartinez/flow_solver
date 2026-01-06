import logging
import torch
import torch.nn as nn
from dataclasses import dataclass
from networks.base_model import BaseModel
from hyperparameters import Hyperparameters


@dataclass
class TrainResult:
    training_losses: list[float]
    validation_losses: list[float]


class Trainer():
    def __init__(self,
                 optimizer: torch.optim.Optimizer,
                 criterion: torch.nn.modules.loss._Loss,
                 training_dataloader: torch.utils.data.DataLoader,
                 validation_dataloader: torch.utils.data.DataLoader = None,
                 ):
        self.optimizer = optimizer
        self.criterion = criterion
        self.training_dataloader = training_dataloader
        self.validation_dataloader = validation_dataloader
        self.logger = logging.getLogger(__name__)

    def train(self, network: BaseModel, epochs: int) -> TrainResult:
        self.logger.info(f"Started training network {network.get_name()} with params: \n{network.params}")

        training_losses: list[float] = []
        validation_losses: list[float] = []

        max_validation_loss = 1_000_000
        for epoch in range(epochs):
            network.train(True)
            training_loss = self._train_one_epoch(network)
            training_losses.append(training_loss)

            network.eval()
            validation_loss = self._validate_one_epoch(network)
            validation_losses.append(validation_loss)

            if (epoch + 1) % (min(max(epochs * 0.20, 1), 200)) == 0:
                self.logger.info(
                    f"Epoch {(epoch + 1):4d} / {epochs:4d}, Train Loss: {training_loss:.4e}, Validation Loss: {validation_loss:.4e}")

            # We only want to export the network with the least loss
            if (validation_loss < max_validation_loss):
                self.logger.info(f"Exporting model {network.get_name()} at epoch {(epoch + 1):4d} with validation loss {validation_loss:.4e}")
                network.export_model()
                max_validation_loss = validation_loss

        self.logger.info(f"Finished training network {network.get_name()}")
        return TrainResult(training_losses, validation_losses)

    def _train_one_epoch(self, network: nn.Module):
        running_loss = 0
        for input, target in self.training_dataloader:
            # Ensure data is float32
            input = input.float()
            target = target.float()

            self.optimizer.zero_grad()

            outputs = network(input)
            loss = self.criterion(outputs, target)
            running_loss += loss.item()
            loss.backward()

            self.optimizer.step()

        return running_loss / len(self.training_dataloader)

    def _validate_one_epoch(self, network: nn.Module):
        running_loss = 0
        with torch.no_grad():
            for input, target in self.validation_dataloader:
                input = input.float()
                target = target.float()

                output = network(input)
                loss = self.criterion(output, target)
                running_loss += loss.item()

        return running_loss / len(self.validation_dataloader)


def train_all_models(training_dataloader: torch.utils.data.DataLoader,
                     validation_dataloader: torch.utils.data.DataLoader,
                     model_export_path: str) -> dict[str, TrainResult]:

    results = {}
    # Loop over all defined networks and train them then save them
    for network in BaseModel.subclasses:
        net = network(model_export_path)
        params = net.params
        # TODO: The optimizer and criterion should be changable
        optimizer = torch.optim.Adam(net.parameters(), lr=params.learning_rate)
        criterion = torch.nn.MSELoss()
        trainer = Trainer(optimizer,
                          criterion,
                          training_dataloader,
                          validation_dataloader)

        result = trainer.train(net, params.epochs)
        results[net.get_name()] = result

    return results
