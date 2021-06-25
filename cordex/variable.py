# -*- coding: utf-8 -*-
# flake8: noqa
"""Variable module

This module defines variables for the CORDEX data request.

"""

from .tables import variables as TABLES


def table(name):
    """Top level function that returns a CORDEX data request table.

    Args:
      name (str): name of the CORDEX table.

    Returns:
      table (DataFrame): Cordex table.

    """
    return TABLES[name]


def tables():
    """Top level function that returns a list of available CORDEX data request tables.

    Returns:
      names (list): list of available CORDEX data request tables.

    """
    return list(TABLES.keys())


class Variable:
    """The :class:`Variable` holds data and meta information of a Cordex Variable.


    **Attributes:**
        *variable_id:*
            cordex variable id.
        *project_id:*
            cordex project id.

    """

    def __init__(self, variable_id, project_id="cmip5"):
        self.variable_id = variable_id
        self.project_id = project_id
        self._from_table()

    def __str__(self):
        return str(self.series)

    def __repr__(self):
        return str(self.series)

    def __getattr__(self, attr):
        if attr in self.series:
            return self.series[attr]
        else:
            raise AttributeError

    def _from_table(self):
        self.series = table(self.project_id).loc[self.variable_id]


def variables(project_id="cmip5"):
    """Top level function that returns all CORDEX variables.

    Args:
      project_id (str): project_id.

    Returns:
      variables (dict): dictionary of Cordex variables.

    """
    return {id: Variable(id, project_id) for id in table(project_id).index}


def variable(variable_id, project_id="cmip5"):
    """Top level function that returns a CORDEX variable.

    Args:
      variable_id (str): name of the variable.

    Returns:
      variable (:class:`_Variable`): Cordex variable.

    """
    return Variable(variable_id, project_id)
