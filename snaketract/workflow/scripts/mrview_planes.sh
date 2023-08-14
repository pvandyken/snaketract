#!/bin/bash
set -euo pipefail
export IFS=$'\n'$'\t'' '

print_help () {
    cat <<EOF
USAGE: mrview_capture img [slice]

Outputs arguments that can be passed to mrview to capture a series of images
across the axial, coronal, and sagital slices

Args
    img     the image to be captured. Used to calculate slice dimensions
    slice   Optional argument specifying the range of captures to be made.
            Follows pythonic START:STOP:STEP syntax. Defaults to ::10
EOF
}

function join_by {
  local d=${1-} f=${2-}
  if shift 2; then
    printf %s "$f" "${@/#/$d}"
  fi
}

params=()
while [[ -n "${1:-}" && ! "$1" == "--" ]]; do
    case "${1:-help}" in
        --help | -h)
            print_help
            exit 0
            ;;
        -* | --*)
            echo "Error: Unsupported flag $1" >&2
            exit 1
            ;;
        * )
            params["${#params[@]}"]="$1"
            ;;
    esac
    shift
done

if [[ "${#params[@]}" -gt 2 ]]; then
    echo "Error: unrecognized args '${params[@]:2}'" >&2
    exit 1
fi

if [[ "${#params[@]}" -eq 0 ]]; then
    echo "Error: must provide path to image" >&2
    exit 1
fi


# parse and validate image
img="${params[0]}"
if [[ ! -e $img ]]; then
    echo "Error: Image does not exist: '$img'" >&2
    exit 1
fi


# parse and validate slice
eval $( echo "${params[1]:-"::"}" | awk -F':' '
    {print "start="$1" stop="$2" step="$3" err="$4}
')
start="${start:-0}"
step="${step:-10}"

if [[ -n $err ]]; then
    echo 'Error: slice must not have more than 3 numbers (START:STOP:SLICE)' >&2
    exit 1
fi

is_number () {
    [[ "$1" =~ ^-?[0-9]+$ ]]
}

for ix in "$start" "${stop:-0}" "${step#"~"}"; do
    is_number "$ix" || (
        >&2 echo "Error: '$ix' in slice '${params[1]:-"::"}' is not a number" &&
        exit 1
    )
done

if [[ -n "$stop" && "$start" -gt "$stop" ]]; then
    echo "Error: Start value must not be greater than stop: '$start > $stop'" >&2
    exit 1
fi


planes=(sagittal coronal axial)

shape=($(mrinfo -size $img))

i_mid=$(( shape[0] / 2 ))
j_mid=$(( shape[1] / 2 ))
k_mid=$(( shape[2] / 2 ))


for plane_i in 0 1 2; do
    step_="${step#"~"}"
    offset=0
    if [ "${step:0:1}" = "~" ]; then
        if [ "$step_" -ge "${shape[$plane_i]}" ]; then
            >&2 echo "$step_ captures can not be grabbed from dimension '${planes[$plane_i]}' with length ${shape[$plane_i]}"
            exit 1
        fi
        step_=$(( shape[$plane_i] / step_ ))
        if [ "$step_" -gt 1 ]; then
            # offset start and stop by half a step to get midway through the
            # image
            offset=$(( step_ / 2 ))
        fi
    fi
    printf -- '-plane %s -capture.prefix %s ' "$plane_i" "${planes[$plane_i]}"
    ceil="${shape[$plane_i]}"
    stop_="${stop:-$ceil}"
    stop_=$(( stop_ + offset ))
    if [[ "$stop_" -gt "$ceil" ]]; then
        stop_="$ceil"
    fi
    i=$(( start + offset ))
    while [[ "$i" -lt "$stop_" ]]; do
        coord=($i_mid $j_mid $k_mid)
        coord[$plane_i]="$i"
        coord_str=$(join_by , ${coord[@]})
        printf -- '-voxel %s -capture.grab ' "$coord_str"
        i=$(( i + step_ ))
    done
done
printf "\n"