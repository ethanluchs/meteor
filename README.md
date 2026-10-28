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

## Dashboard

    detector/   single-file detector (meteor_detect.c), Pi-only build
    ingest/     CSV -> SQLite ingester (Python 3.11 stdlib)
    dashboard/  Next.js 15 app, built here, served from the Pi
    analysis/   offline notebooks/scripts
    data/       CSVs and SQLite files (gitignored)

Develop on the laptop, run on the Pi. Node 22 on both (`.nvmrc`); the Pi
needs 22.13+ so `node:sqlite` works without `--experimental-sqlite`. No
native npm modules (e.g. better-sqlite3): the build is x86_64, the Pi is
aarch64.

### Local development

Pull a snapshot of the Pi's database for real data:

    scp <pi>:<path>/data/dash.db data/dash.db
    cd dashboard
    cp .env.example .env.local
    npm ci
    npm run dev

### Deploy

Build on the laptop. Standalone output leaves out `static/` and `public/`,
so copy them in:

    cd dashboard
    npm run build
    cp -r .next/static .next/standalone/.next/static
    cp -r public .next/standalone/public
    rsync -a --delete .next/standalone/ <pi>:<path>/dashboard/

On the Pi:

    cd <path>/dashboard
    METEOR_DB=<path>/data/dash.db HOSTNAME=0.0.0.0 PORT=3000 node server.js

The standalone bundle already includes its `node_modules`. If you ever
install dependencies on the Pi, use `npm ci`, never `npm install`, so the
lockfile isn't rewritten on a different architecture.
