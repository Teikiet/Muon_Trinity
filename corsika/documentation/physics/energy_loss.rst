Energy Loss in Matter
=====================

This page collects the references used to describe charged-particle energy loss
in matter. The Bethe formula is the baseline, with density, shell, Barkas, and
radiative corrections applied where appropriate.

Bethe Formula and Density Correction
------------------------------------

The 2018 Review of Particle Physics describes the Bethe formula for energy loss
in matter, calculated to order $Z^2$ by quantum methods:

* `Review of Particle Physics: Passage of particles through matter <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/9ee24556574ec08f4088592981bda92e/rpp2018-rev-passage-particles-matter.pdf>`_

The following report provides a related overview:

* `Ziegler, J. Appl. Phys. / Rev. Appl. Phys. 85, 1249-1272 (1999) <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/4a1d125d3db3798fc9718f57c3683247/SRIM-High_Velocity.pdf>`_
* `SRIM <http://www.srim.org/index.htm>`_

Sternheimer's approximation for high-energy density corrections is described in
the Review of Particle Physics. The associated data tables are available in:

* `Sternheimer, Atomic Data and Nuclear Data Tables 30, 261-271 (1984) <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/82e2bf820f20ebfa982035c3eda960f8/1-s2.0-0092640X84900020-main.pdf>`_
* `Groom et al., Atomic Data and Nuclear Data Tables 78, 183-356 (2001) <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/647767824d8c3a2d539b73591d1c6748/1-s2.0-S0092640X01908617-main.pdf>`_
* `PDG atomic and nuclear properties tables <http://pdg.lbl.gov/2018/AtomicNuclearProperties/index.html>`_

Shell and Barkas Corrections
----------------------------

The original Barkas reference provides a simple implementation-oriented
treatment of shell corrections:

* `Barkas reference <https://www.nap.edu/read/20066/chapter/8>`_

For heavy targets with $Z > 56$, a more modern shell-correction treatment is:

* `Bichsel, Physical Review A 46, 5761 (1992) <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/dcdaf0c5f26e0a6df81770028faa36af/PhysRevA.46.5761.pdf>`_

Barkas corrections to the Bethe formula also account for differences between
positive and negative projectiles. An ICRU update is available here:

* `ICRU update (2014) <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/96cdd093e48b84427fc78d0630877434/newstar.pdf>`_

High-Energy Muon Energy Loss
----------------------------

The following work describes improvements to GEANT3 muon $dE/dx$, also used in
CORSIKA 7. Its radiative-loss effects matter primarily above the multi-TeV
range:

* `GEANT3 muon dE/dx improvements <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/0a07e54410651947544636fac5b772c2/1-s2.0-S016890020001007X-main.pdf>`_