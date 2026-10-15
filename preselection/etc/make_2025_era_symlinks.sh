#!/bin/bash
#
# make_2025_era_symlinks.sh
#
# CMS produced no 2025 MC: the Summer24 campaign has to serve both the 2024 and
# 2025 data eras. The analysis therefore needs to process the same skim files
# twice, once per era, each with that era's own calibration and lumi.
#
# RDataFrame identifies a sample by the *path string* given in its input spec,
# so two samples listing the same file collapse into a single identity -- the
# second sample silently vanishes and its events are attributed to the first.
# Giving the 2025 pass its own path avoids that entirely.
#
# This script creates one symlink per matching Run 3 MC dataset directory:
#
#     <dataset>Summer24for2025  ->  <dataset>
#
# The links are RELATIVE, so they resolve identically through the local mount
# (/cmsuf/data/...) and through the xrootd namespace (/store/...). Nothing is
# copied; each link costs a few bytes.
#
# Datasets are discovered by scanning the skim directories, not from a hardcoded
# list, so a newly skimmed sample is picked up without editing this file. Note
# that a dataset entry may itself be a symlink into another skim version (v38,
# v41, ...) rather than a real directory -- those are followed and treated like
# any other dataset, which is why the scan tests "resolves to a directory"
# rather than "is a directory".
#
# Usage:
#     ./make_2025_era_symlinks.sh                  # dry run: show what would happen
#     ./make_2025_era_symlinks.sh --apply          # actually create the links
#     ./make_2025_era_symlinks.sh --undo --apply   # remove them again
#     ./make_2025_era_symlinks.sh --all            # every Run 3 MC dataset, not just the default set
#     ./make_2025_era_symlinks.sh --match REGEX    # restrict to datasets matching REGEX
#
# Safe to re-run: correct links are left alone, and only symlinks whose name
# ends in the suffix are ever created or removed. Real directories are never
# modified.
#
set -euo pipefail

############################ configuration ############################

SKIM_BASE="${SKIM_BASE:-/cmsuf/data/store/user/phchang/skim}"
SKIM_VERSION="${SKIM_VERSION:-VBSVVH_skim_v30}"

# Suffix appended to the dataset directory name. This MUST match
# MC_ERA_CLONE_SUFFIX in dataset_names_ref.py.
SUFFIX="${SUFFIX:-Summer24for2025}"

# Which datasets to clone, as an extended regex matched against the directory
# name. The default is the HT-binned Z(->nunu)+jets and W(->lnu)+jets sets.
# Cloning a dataset makes it appear in 2025 productions, so the default is
# deliberately narrow rather than "every MC sample present"; widen it with
# --match, or drop it entirely with --all.
DEFAULT_MATCH='^(Zto2Nu|WtoLNu)-4Jets_Bin-HT-'
MATCH="$DEFAULT_MATCH"

################################ main #################################

MODE="dryrun"
UNDO="no"

while [ $# -gt 0 ]; do
    case "$1" in
        --apply)    MODE="apply" ;;
        --dry-run)  MODE="dryrun" ;;
        --undo)     UNDO="yes" ;;
        --all)      MATCH="." ;;
        --match)    shift; [ $# -gt 0 ] || { echo "ERROR: --match needs a regex" >&2; exit 1; }
                    MATCH="$1" ;;
        -h|--help)  sed -n '2,40p' "$0"; exit 0 ;;
        *)  echo "ERROR: unknown option '$1'" >&2
            echo "Usage: $0 [--apply] [--undo] [--all | --match REGEX]" >&2
            exit 1 ;;
    esac
    shift
done

ROOT_DIR="$SKIM_BASE/$SKIM_VERSION"
if [ ! -d "$ROOT_DIR" ]; then
    echo "ERROR: skim directory not found: $ROOT_DIR" >&2
    exit 1
fi

# Run 3 MC skim sets. Data is deliberately excluded: data already carries its
# own real year and needs no era clone.
SKIM_DIRS=()
for d in "$ROOT_DIR"/Run3_*; do
    [ -d "$d" ] || continue
    name=$(basename "$d")
    case "$name" in
        *_Data_*) continue ;;
    esac
    SKIM_DIRS+=("$name")
done

if [ ${#SKIM_DIRS[@]} -eq 0 ]; then
    echo "ERROR: no Run 3 MC skim sets found under $ROOT_DIR" >&2
    exit 1
fi

echo "Skim root : $ROOT_DIR"
echo "Suffix    : $SUFFIX"
echo "Match     : $MATCH"
echo "Mode      : $([ "$UNDO" = yes ] && echo undo || echo create) / $MODE"
if [ "$MODE" = "dryrun" ]; then
    echo "            (nothing will be changed; re-run with --apply)"
fi
echo

n_created=0; n_ok=0; n_removed=0; n_conflict=0

for skim_dir in "${SKIM_DIRS[@]}"; do
    dir="$ROOT_DIR/$skim_dir"
    d_created=0; d_ok=0; d_removed=0

    # Every entry that is not itself a clone and resolves to a directory. A
    # dataset may be a real directory or a symlink into another skim version,
    # so test with -d (which follows) after excluding the suffix by name.
    while IFS= read -r ds; do
        [ -z "$ds" ] && continue
        case "$ds" in
            *"$SUFFIX") continue ;;
        esac
        [ -d "$dir/$ds" ] || continue
        printf '%s\n' "$ds" | grep -Eq "$MATCH" || continue

        link="$dir/${ds}${SUFFIX}"

        if [ "$UNDO" = "yes" ]; then
            if [ -L "$link" ]; then
                [ "$MODE" = "apply" ] && rm -f "$link"
                d_removed=$((d_removed+1))
            fi
            continue
        fi

        # Never clobber a real file or directory sitting at the link path.
        if [ -e "$link" ] && [ ! -L "$link" ]; then
            echo "  ERROR: refusing to touch non-symlink: $skim_dir/${ds}${SUFFIX}" >&2
            n_conflict=$((n_conflict+1))
            continue
        fi

        # Already correct? leave it alone.
        if [ -L "$link" ] && [ "$(readlink "$link")" = "$ds" ]; then
            d_ok=$((d_ok+1))
            continue
        fi

        [ "$MODE" = "apply" ] && ln -sfn "$ds" "$link"   # relative target, on purpose
        d_created=$((d_created+1))
    done < <(ls -1 "$dir")

    if [ "$UNDO" = "yes" ]; then
        printf "  %-28s remove %4d\n" "$skim_dir" "$d_removed"
        n_removed=$((n_removed+d_removed))
    else
        printf "  %-28s new %4d   already ok %4d\n" "$skim_dir" "$d_created" "$d_ok"
        n_created=$((n_created+d_created))
        n_ok=$((n_ok+d_ok))
    fi
done

echo
if [ "$UNDO" = "yes" ]; then
    if [ "$MODE" = "apply" ]; then echo "Removed $n_removed symlink(s)."
    else echo "Would remove $n_removed symlink(s). Re-run with --apply to make the changes."; fi
else
    if [ "$MODE" = "apply" ]; then echo "Created $n_created symlink(s); $n_ok already correct."
    else echo "Would create $n_created symlink(s); $n_ok already correct."
         echo "Re-run with --apply to make the changes."; fi
fi

if [ "$n_conflict" -gt 0 ]; then
    echo "WARNING: $n_conflict path(s) skipped because a real file/dir is in the way." >&2
    exit 1
fi

if [ "$MODE" = "apply" ]; then
    echo
    actual=$(find "$ROOT_DIR" -mindepth 2 -maxdepth 2 -type l -name "*${SUFFIX}" 2>/dev/null | wc -l)
    echo "Symlinks matching *${SUFFIX} now under $SKIM_VERSION: $actual"
fi
