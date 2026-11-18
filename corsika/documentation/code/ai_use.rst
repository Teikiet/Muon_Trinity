AI-Assisted Development
=======================

Artificial-intelligence tools, including code assistants and large language
models, may assist CORSIKA 8 work. They do not replace the contributor's
technical, scientific, legal, or ethical responsibility for a contribution.
Every submitted change remains the responsibility of the human contributor who
opens the merge request.

Appropriate Uses
----------------

AI tools may be used to support routine work, for example to:

* explain existing code or documentation;
* draft tests, documentation, refactorings, or build commands;
* identify potential defects, missing cases, or review questions; and
* propose implementation alternatives for discussion.

Treat all generated output as an untrusted draft. Read and understand it before
adding it to the repository, and adapt it to the project's coding style,
architecture, and documentation standards.

Verification and Scientific Responsibility
------------------------------------------

Do not rely on an AI tool as evidence that code, a physical model, a numerical
method, a reference, or a scientific statement is correct. The contributor
must:

* verify the behavior against primary sources, established theory, or an
  appropriate domain expert;
* add and run suitable tests, and provide physics validation when the change
  affects a model, simulation result, or release scope;
* check units, numerical precision, boundary conditions, error handling, and
  performance implications where relevant; and
* verify all generated citations, quotations, links, package versions, and
  licensing statements against their original sources.

The :ref:`release-versioning` requirements for physics validation and feature
freeze apply equally to AI-assisted contributions.

Data, Security, and Confidentiality
-----------------------------------

Use only AI services and configurations approved for the data being processed.
Do not submit the following to an external AI service unless the service and
project owners explicitly permit it:

* credentials, tokens, private keys, passwords, or configuration secrets;
* unpublished research data, simulation results, manuscripts, reviewer
  comments, or collaboration material;
* personal data or information subject to contractual, institutional, or legal
  restrictions; or
* third-party source code, data, or documentation whose licence or agreement
  does not permit that disclosure.

Minimize every prompt to the information required for the task. Remove secrets
from logs, configuration files, issue descriptions, and code snippets before
sharing them with any tool.

Copyright, Licensing, and Attribution
-------------------------------------

AI-generated output can reproduce or closely resemble third-party material.
Before submitting generated code, text, data, or documentation, check that it
is compatible with the CORSIKA 8 License and any relevant module
or external dependency licence. Do not submit copied output whose provenance
or licence cannot be established.

Do not use AI-generated citations as a substitute for reading and citing the
original scientific work. Follow the
`usage and collaboration agreement <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/blob/main/usage_and_collaboration_agreement.md>`_
for project and third-party attribution requirements.

Merge Requests and Review
-------------------------

AI-assisted changes follow the same issue, review, CI, and testing requirements
as all other changes. The merge-request author must be able to explain every
part of the submitted change and respond to reviewer questions without relying
on the tool.

Disclose material AI assistance in the merge request when it significantly
contributed to implementation, analysis, scientific text, or documentation.
State the affected files or task and the verification performed. Routine editor
completion, spelling correction, and formatting assistance need not be called
out separately.

Reviewers should evaluate the submitted change on its technical and scientific
merits. They may request additional tests, source verification, attribution, or
human-authored explanation when the origin, correctness, or licence status is
unclear.
