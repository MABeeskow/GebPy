#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		config.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		16.10.2025

#-----------------------------------------------

"""
Module: config.py
Defines configuration parameters for synthetic mineral data generation.
Used by analysis and synthesis modules within gebpy.core.minerals.
"""

# PACKAGES
from typing import Optional

# CODE
class MineralConfiguration:
    """
    Configuration of the mineral data generation.
    - name: mineral name (e.g., 'Olivine')
    - n_datapoints: Number of generated data points (> 0)
    - random_seed: integer seed; if None, defaults to 42
    """

    def __init__(
            self,
            name: str,
            n_datapoints: int,
            random_seed: Optional[int] = None
    ) -> None:
        self.name = name
        self.n_datapoints = n_datapoints
        self.random_seed = 42 if random_seed is None else random_seed
        self._validate()

    def __repr__(self) -> str:
        return (
        f"<MineralConfiguration name={self.name!r}, "
        f"n_datapoints={self.n_datapoints}, "
        f"random_seed={self.random_seed}>"
        )

    def set_name(self, new_name: str) -> None:
        if not new_name:
            raise ValueError("Mineral name must not be empty. Please assign a valid mineral name.")
        self.name = new_name

    def set_number_of_datapoints(self, new_n_datapoints: int) -> None:
        if not isinstance(new_n_datapoints, int):
            raise TypeError("n_datapoints must be an integer.")
        if new_n_datapoints <= 0:
            raise ValueError("The number of data points must be positive.")
        self.n_datapoints = new_n_datapoints

    def set_random_seed(self, new_random_seed: Optional[int]) -> None:
        self.random_seed = 42 if new_random_seed is None else new_random_seed
        if self.random_seed < 0:
            raise ValueError("The random number seed must be >= 0 or None.")

    def _validate(self) -> None:
        if not self.name:
            raise ValueError("Mineral name must not be empty. Please assign a valid mineral name.")
        if not isinstance(self.n_datapoints, int):
            raise TypeError("n_datapoints must be an integer.")
        if self.n_datapoints <= 0:
            raise ValueError("The number of data points must be positive.")
        if self.random_seed < 0:
            raise ValueError("The random number seed must be >= 0 or None.")

# DEFAULT EXAMPLE
DEFAULT_CONFIG = MineralConfiguration(name="Olivine", n_datapoints=10)