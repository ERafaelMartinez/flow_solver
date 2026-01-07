import logging
import torch
from math import sqrt
from dataclasses import dataclass
from networks.base_model import BaseModel
from normalizer import Normalizer


@dataclass
class TestResult:
    loss_rmse_u: float = 0
    loss_rmse_v: float = 0
    loss_mae_u: float = 0
    loss_mae_v: float = 0
    accuracy: float = 0


class Tester():
    def __init__(self,
                 testing_dataloader: torch.utils.data.DataLoader,
                 normalizer: Normalizer,
                 ):
        self.criterion_mse = torch.nn.MSELoss()
        self.criterion_mae = torch.nn.L1Loss()
        self.testing_dataloader = testing_dataloader
        self.normalizer = normalizer
        self.logger = logging.getLogger(__name__)

    def test(self, network: BaseModel) -> TestResult:
        self.logger.info(f"Started testing network {network.get_name()}")

        # Load last parameters of the network
        network.load_model()
        network.eval()

        running_loss_mse_u = 0
        running_loss_mse_v = 0
        running_loss_mae_u = 0
        running_loss_mae_v = 0
        correct = 0
        total = 0
        with torch.no_grad():
            for input, target in self.testing_dataloader:
                input = input.float()
                target = target.float()

                output = network(input)

                _, target_denormalized = self.normalizer.denormalize(
                    None, target)
                _, output_denormalized = self.normalizer.denormalize(
                    None, output)

                target_denormalized_u = target_denormalized[0, 0, :, :]
                target_denormalized_v = target_denormalized[0, 1, :, :]
                output_denormalized_u = output_denormalized[0, 0, :, :]
                output_denormalized_v = output_denormalized[0, 1, :, :]

                # MSE loss per channel
                loss_mse_u = self.criterion_mse(output_denormalized_u,
                                                target_denormalized_u)
                running_loss_mse_u += loss_mse_u.item()

                loss_mse_v = self.criterion_mse(output_denormalized_v,
                                                target_denormalized_v)
                running_loss_mse_v += loss_mse_v.item()

                # MAE loss per channel
                loss_mae_u = self.criterion_mae(output_denormalized_u,
                                                target_denormalized_u)
                running_loss_mae_u += loss_mae_u.item()

                loss_mae_v = self.criterion_mae(output_denormalized_v,
                                                target_denormalized_v)
                running_loss_mae_v += loss_mae_v.item()

                predictions = output_denormalized.argmax(dim=1)
                correct += (predictions == target_denormalized).sum().item()
                total += target_denormalized.size(0)

        n_testing_dataset = len(self.testing_dataloader)

        test_loss_rmse_u = sqrt(running_loss_mse_u / n_testing_dataset)
        test_loss_rmse_v = sqrt(running_loss_mse_v / n_testing_dataset)
        test_loss_mae_u = running_loss_mae_u / n_testing_dataset
        test_loss_mae_v = running_loss_mae_v / n_testing_dataset

        test_accuracy = correct / total

        result = TestResult(test_loss_rmse_u,
                            test_loss_rmse_v,
                            test_loss_mae_u,
                            test_loss_mae_v,
                            test_accuracy)

        self.logger.info(f"Finished testing network {network.get_name()} with result\n{result}")

        return result


def test_all_models(testing_dataloader: torch.utils.data.DataLoader,
                    model_export_path: str,
                    normalizer: Normalizer) -> dict[str, TestResult]:
    # Loop over all defined networks and test them agains the testing data
    results = {}
    for network in BaseModel.subclasses:
        net = network(model_export_path)
        tester = Tester(testing_dataloader, normalizer)
        result = tester.test(net)
        results[net.get_name()] = result

    return results
