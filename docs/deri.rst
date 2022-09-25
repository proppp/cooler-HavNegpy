HN derivative Module
====================================


``HN_derivative class``
------------------------------------------------------------

.. autoclass:: HavNegpy.HN_derivative
   

``HavNegpy.HN_derivative fit functions``
------------------------------------------------------------


List of fit functions available in the derivative_HN module.

.. automethod:: HavNegpy.HN_derivative.deri_hn

.. automethod:: HavNegpy.HN_derivative.deri_hn_ep

.. automethod:: HavNegpy.HN_derivative.deri_double_hn

.. automethod:: HavNegpy.HN_derivative.ep_s


``HavNegpy.HN_derivative dump methods``
------------------------------------------------------------

Dumps initial guess parameters using the relevant dump methods.


.. automethod:: HavNegpy.HN_derivative.dump_parameters_deri_hn

.. automethod:: HavNegpy.HN_derivative.dump_parameters_deri_double_hn



``HavNegpy.HN_derivative initial view methods``
------------------------------------------------------------

The initial view methods display the initial fit based on the guess parameters supplied.

.. automethod:: HavNegpy.HN_derivative.initial_view_deri_hn

.. automethod:: HavNegpy.HN_derivative.initial_view_deri_hn_ep

.. automethod:: HavNegpy.HN_derivative.initial_view_deri_double_hn


``HavNegpy.HN_derivaitve methods for fitting and saving the fit results``
---------------------------------------------------------------------

The fit method performs the final fit and the final parameters can be saved into a file using save fit method.


.. automethod:: HavNegpy.HN_derivative.create_analysis_file

.. automethod:: HavNegpy.HN_derivative.select_range

.. automethod:: HavNegpy.HN_derivative.fit

.. automethod:: HavNegpy.HN_derivative.save_fit_deri_hn

.. automethod:: HavNegpy.HN_derivative.save_fit_deri_hn_ep

.. automethod:: HavNegpy.HN_derivative.save_fit_deri_double_HN


