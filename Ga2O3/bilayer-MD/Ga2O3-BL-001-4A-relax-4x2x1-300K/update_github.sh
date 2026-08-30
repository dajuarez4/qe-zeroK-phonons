#!/bin/bash

set -euo pipefail

PROJECT="/Users/dajuarez4/Documents/qe-zeroK-phonons/Ga2O3/bilayer-MD/Ga2O3-BL-001-4A-relax-4x2x1-300K"
REMOTE="${GIT_REMOTE:-origin}"
MAX_FILE_MB="${MAX_FILE_MB:-25}"

cd "$PROJECT"

if [[ ! "$MAX_FILE_MB" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MAX_FILE_MB debe ser un número entero positivo."
    exit 1
fi

REPO_ROOT="$(git rev-parse --show-toplevel)"
PROJECT_REL="$(git rev-parse --show-prefix)"
PROJECT_REL="${PROJECT_REL%/}"
BRANCH="$(git branch --show-current)"

if [[ -z "$BRANCH" ]]; then
    echo "ERROR: Git está en detached HEAD; selecciona una rama antes de continuar."
    exit 1
fi

if ! git remote get-url "$REMOTE" >/dev/null 2>&1; then
    echo "ERROR: no existe el remote '$REMOTE'."
    exit 1
fi

# No mezclar cambios que el usuario ya hubiera dejado preparados.
if ! git diff --cached --quiet; then
    echo "ERROR: ya existen cambios staged en Git."
    echo "Haz commit o unstage de esos cambios antes de usar este script."
    git status --short
    exit 1
fi

echo "Repositorio: $REPO_ROOT"
echo "Proyecto:    $PROJECT_REL"
echo "Rama:        $BRANCH"
echo "Límite:      ${MAX_FILE_MB} MB por archivo"
echo

# Agregar únicamente este proyecto, no otros cálculos del repositorio.
git -C "$REPO_ROOT" add -A -- "$PROJECT_REL"

echo "Archivos de simulación excluidos:"
excluded_count=0
while IFS= read -r -d '' path; do
    local_path="${path#"$PROJECT_REL"/}"
    case "$local_path" in
        Ga2O3-BL-001-4x2x1.md.out|\
        TDEP_300K/infile.forces|\
        TDEP_300K/infile.positions|\
        TDEP_300K/outfile.dispersion_relations.hdf5|\
        HELD_300K/held_step_coefficients.csv|\
        HELD_300K/held_step_dashboard_cache.npz|\
        HELD_300K/Ga2O3-BL-001-4A-HELD-step-dashboard.gif|\
        *.hdf5|*.npz)
            git -C "$REPO_ROOT" restore --staged -- "$path"
            echo "  - $local_path"
            excluded_count=$((excluded_count + 1))
            ;;
    esac
done < <(git -C "$REPO_ROOT" diff --cached --name-only -z)

# Segunda barrera: excluir cualquier blob staged que exceda el límite.
max_bytes=$((MAX_FILE_MB * 1024 * 1024))
while IFS= read -r -d '' path; do
    size="$(git -C "$REPO_ROOT" cat-file -s ":$path")"
    if (( size > max_bytes )); then
        git -C "$REPO_ROOT" restore --staged -- "$path"
        size_mb=$(( (size + 1024 * 1024 - 1) / (1024 * 1024) ))
        echo "  - ${path#"$PROJECT_REL"/} (${size_mb} MB; supera el límite)"
        excluded_count=$((excluded_count + 1))
    fi
done < <(
    git -C "$REPO_ROOT" diff --cached \
        --name-only --diff-filter=ACMR -z -- "$PROJECT_REL"
)

if (( excluded_count == 0 )); then
    echo "  (ninguno)"
fi

if git -C "$REPO_ROOT" diff --cached --quiet -- "$PROJECT_REL"; then
    echo
    echo "No hay cambios ligeros para subir."
    exit 0
fi

echo
echo "Cambios que SÍ se incluirán en el commit:"
git -C "$REPO_ROOT" diff --cached --stat -- "$PROJECT_REL"

if ! git -C "$REPO_ROOT" diff --cached --check -- "$PROJECT_REL"; then
    echo "ERROR: Git encontró problemas en los cambios staged."
    exit 1
fi

echo
read -r -p "¿Crear el commit y hacer push a $REMOTE/$BRANCH? [y/N] " answer
case "$answer" in
    y|Y|yes|YES|s|S|si|SI|sí|SÍ) ;;
    *)
        git -C "$REPO_ROOT" restore --staged -- "$PROJECT_REL"
        echo "Cancelado. Tus archivos locales no fueron modificados."
        exit 0
        ;;
esac

if (( $# > 0 )); then
    commit_message="$*"
else
    commit_message="Update TDEP results $(date '+%Y-%m-%d %H:%M')"
fi

# El índice comenzó vacío y arriba sólo agregamos este proyecto; por eso
# hacemos commit del índice tal como fue revisado, sin volver a leer archivos.
git -C "$REPO_ROOT" commit -m "$commit_message"
git -C "$REPO_ROOT" pull --rebase --autostash "$REMOTE" "$BRANCH"
git -C "$REPO_ROOT" push "$REMOTE" "$BRANCH"

echo
echo "GitHub actualizado correctamente en $REMOTE/$BRANCH."
