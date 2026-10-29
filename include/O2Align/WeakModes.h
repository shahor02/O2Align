// Diagnostics of the global weak (deformation) modes of an alignment result
#ifndef O2ALIGN_WEAKMODES_H
#define O2ALIGN_WEAKMODES_H

#include <array>
#include <map>
#include <string>
#include <cstdint>
#include <vector>

namespace o2::alignrs
{
class Volume;
class Detector;

/// Displacement of a sensor produced by a set of fitted (or injected) rigid-body values.
struct SensorDisplacement {
  const Volume* top{};        ///< top volume of the detector branch the sensor belongs to
  const Volume* vol{};        ///< the sensor
  int layer{-1};              ///< layer parsed from the symbolic name, -1 if none
  std::array<double, 3> x{};  ///< position of the sensor centre, cm
  std::array<double, 3> d{};  ///< displacement of the centre, cm
  std::array<double, 3> w{};  ///< rotation vector, rad
};

/// Compose the rigid-body values (raw label -> value, absent = 0) down the hierarchy, like
/// Detector::MP2AlignParams does, into the displacement of every sensor, i.e. of every volume which is a
/// leaf and neither virtual nor pseudo. An ITS3 half barrel (a volume whose children are all pseudo tiles)
/// is a sensor as well: the displacement of its local origin and its rotation are given, the non-rigid
/// calibration DOFs are not included. Independent of the hierarchy level which carries the values, which
/// makes it the quantity to compare an alignment result with the injected misalignment of a closure test.
std::vector<SensorDisplacement> computeSensorDisplacements(Volume* root, const std::map<uint32_t, double>& labelToValue);

/// Diagnose the global deformation modes of Volume::GlobalMode in the fitted rigid-body corrections.
///
/// A weak mode is a coherent deformation of a whole detector branch which changes the sensor positions
/// but (nearly) not the track chi2: a track from the vertex absorbs it into its own parameters (1/pt,
/// phi, tgLambda, d0, z0). The hierarchical mean constraint removes only net rigid motions of a parent's
/// children, so these modes are not constrained by it. The modes below are the ones implemented, as
/// displacement fields in cylindrical coordinates (r, phi, z), t = z/r, eps the amplitude. They can be
/// constrained in the fit by the "globalModes" clause of the DOF config (Volume::writeGlobalModeConstraints,
/// doc/DOFConfig_rules.md); each is diagnosed here by one column per equation, named as in the clause.
///
///  mode (equations)                  field                                   what it is
///  radialScaling (1)                 dr = eps*r                              uniform transverse dilation; absorbed by 1/pt (B!=0)
///  zScaling (1)                      dz = eps*z                              longitudinal scale; absorbed by tgLambda
///  curl (1)                          r*dphi = eps*r                          rigid rotation about the beam (also a rigid DOF)
///  twist (1)                         r*dphi = eps*z*r                        torsion: the two ends rotated oppositely
///  elliptical (2)                    dr = eps*cos2phi, eps*sin2phi           layers squashed into ellipses; r-independent amplitude
///  bowing (1)                        dr = eps*z^2                            radial sagitta along the beam
///  radialHarmonic<l> (2 per l)       dr = eps*r*cos(l*phi), eps*r*sin(l*phi)         l = 1..4; l=0 is radialScaling
///  tangentialHarmonic<l> (2 per l)   r*dphi = eps*r*cos(l*phi), eps*r*sin(l*phi)     l = 1..4; l=0 is curl
///  globalTilt (2: _x, _y)            rigid rotation about the global X / Y   transverse partners of curl
///  radialPolar_r<p>t<k> (5)          dr = eps*r^p*t^k                        (p,k) = (0,0),(0,1),(0,2),(1,1),(1,2)
///  tangentialPolar_r<p>t<k> (5)      r*dphi = eps*r^p*t^k                    same (p,k)
///  longitudinalPolar_r<p>t<k> (4)    dz = eps*r^p*t^k                        (0,1),(0,2),(1,0),(1,2)
///
/// The polar families are the fields depending on the polar angle: a track from the vertex is at z = r*t
/// at radius r, so r*f(t) is absorbed by its phi/tgLambda/curvature (null directions whatever the vertex
/// does), f(t) only by d0/z0 (broken by the vertex). Their duplicates of other modes are left out. The
/// modes with a handle only on DOFs distributed in z (zScaling, bowing, polar families) are not
/// observable from stave-level DOFs, all centred at z = 0. In total Volume::NGlobalEq = 42 columns. The
/// globalShift columns duplicate the rigid shift of the fit, like curl and globalTilt its rotation, so
/// their B amplitudes carry no information; they exist for the "groupModes" constraints.
///
/// Not diagnosed here, but available as a constraint: the "rayModes" clause, which removes the common
/// angular displacement r*F(phi, t) of the sensors of every (phi, t) cell, i.e. the whole function space
/// of which the harmonic and polar modes above are low-order members.
///
/// For every detector branch (child of the root) two complementary numbers are given per mode:
///  - metric A (what Millepede's mode constraint sees): sum(c_i * p_i) / sum(|c_i|) over the free
///    rigid-body DOFs of the branch, c_i being the coefficient the DOF would have in the mode equation
///    written by Volume::writeGlobalModeConstraints, p_i its fitted value. In micrometers (rotations
///    enter with the lever arm of the mode field).
///  - metric B (what the geometry actually did): the corrections of the whole hierarchy are composed
///    into the displacement of the centre of every sensor, which is fitted by a global rigid shift +
///    rotation plus the mode fields. The amplitudes are given as the displacement the mode produces at
///    the largest lever arm of the branch, with the fraction of the displacement power they explain and
///    the mode-to-mode correlations. Independent of which hierarchy level carries the correction.
/// Writes <outBase>.json (numbers), <outBase>.root (per-sensor displacements and fit results, see
/// macro/drawWeakModes.C) and logs a summary table. The fitted values are absolute: an injected
/// misalignment (misAlgJson) is not subtracted.
void diagnoseWeakModes(Volume* root, const std::map<uint32_t, double>& labelToValue, const std::string& outBase);

/// Diagnose the global deformation modes already present in the initial alignment itself, i.e. the
/// CCDB o2::detectors::AlignParam objects baked into the loaded geometry (Detector::MP2AlignParams'
/// "step 2"), relative to the ideal geometry. No millepede fit is involved: metric A (which needs a
/// fit's free-DOF values) is not computed, only metric B (the sensor displacement decomposition).
/// Same outputs as diagnoseWeakModes, minus the metricA_um/nTermsA fields.
void diagnoseInitialAlignment(Volume* root, const std::vector<Detector*>& detectors, const std::string& outBase);

} // namespace o2::alignrs

#endif
