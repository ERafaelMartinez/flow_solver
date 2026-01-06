import csv
import os
import logging
from dataclasses import asdict, fields
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from trainer import TrainResult
from tester import TestResult


class Reporter():
    def __init__(self, reports_export_path: str):
        self.reports_export_path = reports_export_path
        self.logger = logging.getLogger(__name__)

    def create_training_report(self, train_results: dict[str, TrainResult]):
        self.logger.info("Started Creating training reports")
        with PdfPages(os.path.join(self.reports_export_path, "training_report.pdf")) as pdf:
            for network_name, train_result in train_results.items():
                self.logger.debug(
                    f"Creating training report for {network_name}")
                fig, ax = plt.subplots()
                self._plot_training_result(ax, network_name, train_result)
                pdf.savefig(fig)
        self.logger.info("Finished Creating training reports")

    def create_testing_report(self, test_results: dict[str, TestResult]):
        self.logger.info("Started Creating testing reports")
        self.logger.debug("Exporting raw testing data")
        test_results_to_csv = [
            {'name': k, **asdict(v)} for (k, v) in test_results.items()
        ]
        with open(os.path.join(self.reports_export_path, 'test_results_raw.csv'), 'w', newline='') as csvfile:
            writer = csv.DictWriter(
                csvfile, fieldnames=test_results_to_csv[0].keys())
            writer.writeheader()
            writer.writerows(test_results_to_csv)

        # We want to generate a bar chart per metric in the TestResult dataclass
        # The goal is to automatically do this here so that we developers only
        # need to worry about the implementation of the tester class 
        # and not about the graphs
        self.logger.debug("Exporting pdf testing report")
        with PdfPages(os.path.join(self.reports_export_path, "testing_report.pdf")) as pdf:
            for field in fields(TestResult):
                field_results = {k: v.__getattribute__(field.name) for k, v in test_results.items()}
                fig, ax = plt.subplots()
                self._plot_testing_results(ax, field_results, field.name)
                pdf.savefig(fig)

        self.logger.info("Finished Creating testing reports")

    def _plot_testing_results(self,
                              ax: plt.Axes,
                              results: dict[str, float],
                              yLabel: str):
        ax.bar(results.keys(), results.values(), label=f"Test {yLabel}")
        ax.legend()
        ax.set_xlabel("Networks")
        ax.set_ylabel(yLabel)
        ax.set_title(f"Testing {yLabel}")

    def _plot_training_result(self,
                              ax: plt.Axes,
                              model_name: str,
                              train_result: TrainResult):
        ax.semilogy(train_result.training_losses, label="Train loss")
        ax.semilogy(train_result.validation_losses, label="Validation loss")
        ax.legend()
        ax.set_xlabel("Epoch")
        ax.set_ylabel("Loss (log scale)")
        ax.set_title(f"Training and validation loss ({model_name})")
