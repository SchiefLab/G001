"""
This type stub file was generated by pyright.
"""

from ._core import VectorPlotter
from ._decorators import _deprecate_positional_args

__all__ = ["relplot", "scatterplot", "lineplot"]
_relational_narrative = ...
_relational_docs = ...
_param_docs = ...
class _RelationalPlotter(VectorPlotter):
    wide_structure = ...
    sort = ...
    def add_legend_data(self, ax):
        """Add labeled artists to represent the different plot semantics."""
        ...
    


class _LinePlotter(_RelationalPlotter):
    _legend_attributes = ...
    _legend_func = ...
    def __init__(self, *, data=..., variables=..., estimator=..., ci=..., n_boot=..., seed=..., sort=..., err_style=..., err_kws=..., legend=...) -> None:
        ...
    
    def aggregate(self, vals, grouper, units=...): # -> tuple[Unknown, Unknown, None] | tuple[Unknown, Unknown, DataFrame | Series[Unknown] | Unknown | None]:
        """Compute an estimate and confidence interval using grouper."""
        ...
    
    def plot(self, ax, kws): # -> None:
        """Draw the plot onto an axes, passing matplotlib kwargs."""
        ...
    


class _ScatterPlotter(_RelationalPlotter):
    _legend_attributes = ...
    _legend_func = ...
    def __init__(self, *, data=..., variables=..., x_bins=..., y_bins=..., estimator=..., ci=..., n_boot=..., alpha=..., x_jitter=..., y_jitter=..., legend=...) -> None:
        ...
    
    def plot(self, ax, kws): # -> None:
        ...
    


@_deprecate_positional_args
def lineplot(*, x=..., y=..., hue=..., size=..., style=..., data=..., palette=..., hue_order=..., hue_norm=..., sizes=..., size_order=..., size_norm=..., dashes=..., markers=..., style_order=..., units=..., estimator=..., ci=..., n_boot=..., seed=..., sort=..., err_style=..., err_kws=..., legend=..., ax=..., **kwargs):
    ...

@_deprecate_positional_args
def scatterplot(*, x=..., y=..., hue=..., style=..., size=..., data=..., palette=..., hue_order=..., hue_norm=..., sizes=..., size_order=..., size_norm=..., markers=..., style_order=..., x_bins=..., y_bins=..., units=..., estimator=..., ci=..., n_boot=..., alpha=..., x_jitter=..., y_jitter=..., legend=..., ax=..., **kwargs):
    ...

@_deprecate_positional_args
def relplot(*, x=..., y=..., hue=..., size=..., style=..., data=..., row=..., col=..., col_wrap=..., row_order=..., col_order=..., palette=..., hue_order=..., hue_norm=..., sizes=..., size_order=..., size_norm=..., markers=..., dashes=..., style_order=..., legend=..., kind=..., height=..., aspect=..., facet_kws=..., units=..., **kwargs): # -> FacetGrid:
    ...

