How to use git
==============

The repository of MulSKIPS is in GitHub_. For each developer, it is safe to clone its own repository.

.. _GitHub: https://github.com/MulSKIPS/MulSKIPS

Each developer has its own devel branch which has to be a copy of the MulSKIPS branch (called **main**). 
Only the command ``git pull main`` will be used to update your **devel** branch.
It is highly recommended not to develop in your **devel** branch!

For each set of modifications, start a new branch with ``git checkout -B newbranch`` from the **devel** one.
Then you do your modifications, commits and so on.

Before the merge request, you have to update your **newbranch** branch from the **devel** branch 
(which will be updated from the main one) with ``git rebase devel`` and solve the conflicts if necessary.

Then go to GitHub and do a merge-request. Finally do not change your **newbranch** branch because the merge request will be
affected. In function of the remarks of the reviewer, you change you **newbranch** branch, do a ``git pull origin newbranch`` to propagate the commits. 

Development of a new functionality
==================================

The first question to ask yourself is the *generality* of the
functionality you are going to implement.
The spirit is to work at the lowest possible level for a given task.
The idea is to make available the functionality also to other tasks of MulSKIPS.

Developers are supposed to document as much as possible all introduced modules, functions and
subroutines together with their input/output variables by means of comment lines in the source files.
The same holds for all jupyter notebook developed.

Create a test for the functionality
-----------------------------------

Once the developed functionality is well tested and integrated in the MulSKIPS suite,
it is fundamental to fix such development work with a regression test.
All regression tests are handled and collected here_.

.. _here: https://github.com/MulSKIPS/MulSKIPS/regression-tests

Code developements and the associated functionalities will be lost without a proper regression test
which runs the code and touches all the developed source files.

Make a tutorial which demonstrates the functionality
----------------------------------------------------

Once the developed functionality lies within MulSKIPS, it is important to build a tutorial which
explains how to make advantage of the new functionality.
As for the regression tests, tutorials in MulSKIPS are built with jupyter notebooks.
In principle each tutorial should report as much as possible details on how to run MulSKIPS with the new
functionality, how to analyse and visualize output data.
The wealth of informations at this stage will make the newly implement features successful.
Once an appropriate tutorial notebook has been written, this should be added to the tutorial directory
``MULSKIPS_ROOT/tutorial-mulskips/.``.

How to write documentation in this package
==========================================

The present documentation has been generated with readthedocs_ and the ``sphinx`` package.
The syntax in the text pages is ``ReStructured Text`` format.
See `this page`_ for examples of the syntax, and `this one`__ for a more comprehensive explanation of
the syntax.

.. _readthedocs: https://docs.readthedocs.io/en/stable/index.html

.. __: http://docutils.sourceforge.net/rst.html

.. _this page: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

