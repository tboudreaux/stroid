# Generating Releases

> This document is likely only of interest to developers wishing to build release files. If you wish to use
> stroid, please refer to the main readme instructions.

## Building
Release files can be built on an arm mac using the following commands:

```bash
./release.sh
```

all release files will be placed in the `releases/` folder.

This will construct release files for arm mac, arm linux, and x86_64 linux. These should be
compatible with most systems (linux more recent than 2020 and macOS more recent than 15.0)

For those who wish to build from source see the main readme instructions.

## Packaging
The auto-generated release files are named according to the following convention:
`stroid-<os>-<arch>`, these are all statically linked so they should
just run on system. To distribute these files we want to package them into a single 
archive file along with an installation script which will select the correct binary
for the host system. There is a script `package_release.sh` which will do this automatically.

Once the release files have been generated, run the following command:

```bash
./package_release.sh
```

## Uploading
Once the release files have been packaged, they can be uploaded to the github releases page. This can be done manually through the github web interface, or
automatically using the `gh` CLI tool. To upload using the CLI tool, first ensure you have it installed and authenticated with your
github account. Then run the following command:

```bash
gh release create <tag> ./stroid-release.tar.gz
```

Alternatively, you can upload the files manually through the github web interface by navigating to the releases page of
the repository and clicking on "Draft a new release". From there, you can upload the `stroid-release.tar.gz` file and
publish the release.