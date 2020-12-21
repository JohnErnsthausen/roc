SRS

1.

Table of Symbols: Do you mean sequence or set? Usually {} is for a set.

I would like to use Rudin's notation as I reference his book later.

See [https://en.wikibooks.org/wiki/Real_Analysis/Sequences] for other notations such as ().

Prof Smith gave the notation <> for sequences P.3.

2.

Remove repeated names in Abbreviations and Acronyms

Remove FORTRAN77: You don't normally mention a programming language in the SRS.

3.

Use \cite{} to cite the paper as part of the text with author names and year of publication.

Use \citep{} to put citation in parenthetical

4.

P2 What would a typical title be for the associated course?

Calculus 2
Complex Analysis 1
Real Analysis

5.

P2 Show how ROC fits into a typical context

6.

What you gave are some potential contexts in which ROC could be used. This information could
go into the previous section. What we want in this section is the background needed by a user
of the software.

A calculus student.
A user of MAPLE or MATLAB.
A developer of a Taylor series method. Reference the previous section.

* 7.

Say what the inputs are p5

*8.

Move GS2, CS3, and GS4 to the VnV plan. They are not goals of the software.

GS2: Compute the Radius of Convergence and the Order of Singularity and The Error

*9.

You should have an assumption related to your scope (top-line) p5

*10.

You have reversed the detailed description for the template. The description should
go after the IM in both cases.

*11.

Does the Radius of Convergence have to be positive?

Yes.

*12.

How will ROC find a scale? If this is part of your requirements, then you have to say what this means.

*13.

You already said Functional Requirements R2, and it is not a requirement.

However keep "ROC should execute as fast as the CC software DRDCV" as a nonfunctional requirement.

*14.

R3 is part of VnV Plan.

*15.

R4 should be moved from a Functional Requirement to a non-functional Requirement.

*16.

Yes ROC has functional requirements. Performance requirements are non-functional requirements.
Accuracy requirements are non-functional requirements. Portability and maintainability are
non-functional; requirements.

TODO
Have a look at the rational I added to the ??? project template on non-functional requirements.

