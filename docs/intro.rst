Introduction
============

``HavNegpy`` is an Python package to analyze the dielectric spectrsocopy data. The package is written specifically to analyze the dielectric loss and real part of AC conductivity.
``HavNeg`` is an acronym for Havriliak and Negami function.
 
The package contains three modules and each module can be instantiated as::

              >>> import HavNegpy as h
              >>> hn = h.HN()
              >>> hn_deri = h.HN_derivative()
              >>> cond = h.Conductivity()
	   

All modules contain same methods which includes:
``selecting range of data``, ``dumping initial fit parameters``, ``perform least squares fitting``, ``creating an analysis file to save fit results``, and ``save the fit results``.
Besides, other method include ``initial view of fit parameters``

A clear description is provided in the tutorial.

Note: Some HTML images aren't rendered in the tutorials page. All tutorial notebooks are available at https://github.com/mkolmang/Tutorials_HavNegpy 



