# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.1] - 2026-04-17

### Fixed
- Fixed logging configuration side effects when idat-tools is imported as a library in other projects
  - Replaced `logging.basicConfig()` with `logging.NullHandler()` to prevent forceful DEBUG level configuration
  - Follows Python logging best practices for libraries to avoid interfering with importing projects' logging setup
