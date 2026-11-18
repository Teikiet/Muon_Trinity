/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/modules/writers/ParticleWriterParquet.hpp>
#include <corsika/modules/writers/WriterOff.hpp>
#include <corsika/framework/core/Step.hpp>

#include <optional>

namespace corsika {

  /**
   * @ingroup Modules
   * @{
   *
   * The ObservationPlane writes PDG codes, kinetic energies, locations and momentum unit
   * vectors of particles with respect to the central point of the plane into its output
   * file. By default, the particles are considered "absorbed" afterwards. You can also
   * set the ObservationPlane as non-absorbing. The plane's normal vector defines the
   * "positive side" from which particles are "caught." Only particles crossing the plane
   * from the positive normal side are recorded. The direction vector `x_dir` should be
   * orthogonal to the normal vector and defines the x/y coordinate system for the output
   *
   * The default output format is parquet.
   *
   * **Note/Limitation:** as discussed in
   * https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/issues/397
   * you cannot put two ObservationPlanes exactly on top of each
   * other. Even if one of them is "permeable". You have to put a
   * small gap in between the two plane in such a scenario, or develop
   * another more specialized output class.
   */
  template <typename TTracking, typename TOutput = ParticleWriterParquet>
  class ObservationPlane : public ContinuousProcess<ObservationPlane<TTracking, TOutput>>,
                           public TOutput {

    using TOutput::write;

  public:
    /**
     * Construct a new Observation Plane object.
     *
     * @tparam TArgs
     * @param plane The plane.
     * @param x_dir The x-direction/axis.
     * @param absorbing Flag to make the plane absorbing.
     * @param padding if absorbing, particles will be placed this far past the observation
     * plane to avoid machine-precision issues. Should choose this value to be small based
     * on the medium
     * @param outputArgs arguments that will be given to the particle writer
     */
    template <typename... TArgs>
    ObservationPlane(Plane const& plane, DirectionVector const& x_dir,
                     bool const absorbing = true, LengthType const padding = 1e-6 * 1_m,
                     TArgs&&... outputArgs);

    ~ObservationPlane() {}

    template <typename TParticle>
    ProcessReturn doContinuous(Step<TParticle>&, bool const stepLimit);

    template <typename TParticle, typename TTrajectory>
    LengthType getMaxStepLength(TParticle const&, TTrajectory const& vTrajectory);

    Plane getPlane() const { return plane_; }
    DirectionVector getXAxis() const { return xAxis_; }
    DirectionVector getYAxis() const { return yAxis_; }

    /**
     * Switch the recorded coordinates to global (i.e. the observation point's
     * coordinate system), translated so that `origin` maps to (0, 0, 0).
     *
     * When enabled, the columns `x, y, z` hold the particle's position in global
     * coordinates relative to `origin`, and `nx, ny, nz` hold the direction cosines
     * in that same global system. When never called (the default), the recorded
     * values are plane-local as before, and `z` is zero.
     */
    void setGlobalOutput(Point const& origin);

    YAML::Node getConfig() const;

  private:
    Plane const plane_;
    DirectionVector const xAxis_;
    DirectionVector const yAxis_;
    bool const deleteOnHit_;
    LengthType const padding_;
    /// If set, record global coordinates relative to this origin instead of
    /// plane-local ones.
    std::optional<Point> globalOutputOrigin_;
  };
  //! @}
} // namespace corsika

#include <corsika/detail/modules/ObservationPlane.inl>