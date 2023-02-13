# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.1]

### Fixes
* Geometry: fix method transfer_matrix_2_euler() in the case of R[2, 2] is really near and superior to 1

## [0.3.0]

### New features

* Linkages: Add new type of Linkages: FrictionlessBevelGearLinkage

### Fixes

* Global: update with new version of volmdlr ( copy argument in rotation)
* Global: Change dessia_common import with the new version
* Linkages: Error in GearSetLinkage equation with helical angle
* Dynamic Positions: Function Babylon script for LineSegment3D (remove in volmdlr with the version 0.8.0)

## [0.1.3]
### Changed
- part_frame to part_global_frame
- is_equivalent on linkage

## [0.1.2]
### Added
- Dynamic positions module

### Debug

### Changed
