# meteor

Meteor forward-scatter detection service.

Listens to an RTL-SDR dongle, detects brief signal spikes caused by
meteor ionization trails reflecting a distant radio beacon, logs each
event to SQLite, and POSTs a notification to a local backend.

## Layout

    src/        C sources, one .c/.h pair per concern
    build/      object files (gitignored, arch-specific)
    config/     example config; real config is gitignored

## Building

Requires: gcc, make, pkg-config, librtlsdr-dev, libsqlite3-dev,
libcurl4-openssl-dev

    make            # build
    make clean      # remove build artifacts

## Development

Built on WSL (x86_64) against a file-backed sample source, deployed to
a Raspberry Pi 5 (aarch64) with the real dongle. Source syncs via git;
binaries never do — each machine compiles its own.
