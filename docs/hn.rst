HN Module 
========================

``HN class``
-------------

.. autoclass:: HavNegpy.HN
   

``HavNegpy.HN fit functions``
-------------------------------

List of fit functions available in the HN module.

.. automethod:: HavNegpy.HN.hn

.. automethod:: HavNegpy.HN.hn_cond

.. automethod:: HavNegpy.HN.hn_flank

.. automethod:: HavNegpy.HN.flank_s

.. automethod:: HavNegpy.HN.cond_s

.. automethod:: HavNegpy.HN.double_hn

.. automethod:: HavNegpy.HN.double_hn_cond


``HavNegpy.HN dump methods``
-----------------------------------------

Dumps initial guess parameters using the relevant dump methods.


.. automethod:: HavNegpy.HN.dump_parameters_hn

.. automethod:: HavNegpy.HN.dump_parameters_double_hn

.. automethod:: HavNegpy.HN.dump_parameters_flank


``HavNegpy.HN initial view methods``
----------------------------------------------

The initial view methods display the initial fit based on the guess parameters supplied.

.. automethod:: HavNegpy.HN.initial_view_hn

.. automethod:: HavNegpy.HN.initial_view_hn_cond

.. automethod:: HavNegpy.HN.initial_view_double_hn

.. automethod:: HavNegpy.HN.initial_view_double_hn_cond

.. automethod:: HavNegpy.HN.initial_view_flank


``HavNegpy.HN methods for fitting and saving the fit results``
------------------------------------------------------------

The fit method performs the final fit and the final parameters can be saved into a file using save fit method.

.. automethod:: HavNegpy.HN.create_analysis_file

.. automethod:: HavNegpy.HN.select_range

.. automethod:: HavNegpy.HN.fit

.. automethod:: HavNegpy.HN.save_fit_hn

.. automethod:: HavNegpy.HN.save_fit_hn_cond

.. automethod:: HavNegpy.HN.save_fit_hn_flank

.. automethod:: HavNegpy.HN.save_fit_double_HN


