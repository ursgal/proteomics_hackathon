#!/usr/bin/env python
import pathlib
import sys
import ursgal


def main(mzml_folder, merged_result):
    params = {
        "isotopic_distribution_tolerance": 5,
        "normalize_intensities": True,
        "integrate_peak_areas": False,
        "only_precursor_charge": False,
        "match_between_runs": True,
        "match_between_runs_RT_window": 0.5,  # maybe 1?
        "require_msms_id": False,
        "bayesian_fold_change": True,
        "bayesian_fold_change_control_condition": "H",  # H = Healthy, A = ALS
        "fold_change_cutoff": 0.1,
        "markov_chain_iterations": 3000,
        "markov_chain_burn_in_iterations": 1000,
        "use_shared_peptides": False,
        "random_seed": 200,
    }
    uc = ursgal.UController(verbose=True, params=params, profile="QExactive+",)

    healthy_patient_files = [
        "TN_CSF_062617_02.mzML",
        "TN_CSF_062617_05.mzML",
        "TN_CSF_062617_09.mzML",
        "TN_CSF_062617_12.mzML",
        "TN_CSF_062617_13.mzML",
        "TN_CSF_062617_14.mzML",
        "TN_CSF_062617_16.mzML",
        "TN_CSF_062617_18.mzML",
        "TN_CSF_062617_19.mzML",
        "TN_CSF_062617_22.mzML",
        "TN_CSF_062617_23.mzML",
        "TN_CSF_062617_24.mzML",
        "TN_CSF_062617_25.mzML",
        "TN_CSF_062617_29.mzML",
        "TN_CSF_062617_33.mzML",
        "TN_CSF_062617_34.mzML",
        "TN_CSF_062617_36.mzML",
        "TN_CSF_062617_37.mzML",
        "TN_CSF_062617_38.mzML",
        "TN_CSF_062617_43.mzML",
        "TN_CSF_062617_45.mzML",
        "TN_CSF_062617_46.mzML",
        "TN_CSF_062617_50.mzML",
        "TN_CSF_062617_51.mzML",
        "TN_CSF_062617_52.mzML",
        "TN_CSF_062617_53.mzML",
        "TN_CSF_062617_54.mzML",
        "TN_CSF_062617_57.mzML",
    ]
    als_patient_files = [
        "TN_CSF_062617_03.mzML",
        "TN_CSF_062617_04.mzML",
        "TN_CSF_062617_06.mzML",
        "TN_CSF_062617_07.mzML",
        "TN_CSF_062617_08.mzML",
        "TN_CSF_062617_10.mzML",
        "TN_CSF_062617_11.mzML",
        "TN_CSF_062617_15.mzML",
        "TN_CSF_062617_17.mzML",
        "TN_CSF_062617_20.mzML",
        "TN_CSF_062617_21.mzML",
        "TN_CSF_062617_26.mzML",
        "TN_CSF_062617_27.mzML",
        "TN_CSF_062617_28.mzML",
        "TN_CSF_062617_30.mzML",
        "TN_CSF_062617_31.mzML",
        "TN_CSF_062617_32.mzML",
        "TN_CSF_062617_35.mzML",
        "TN_CSF_062617_39.mzML",
        "TN_CSF_062617_40.mzML",
        "TN_CSF_062617_41.mzML",
        "TN_CSF_062617_42.mzML",
        "TN_CSF_062617_44.mzML",
        "TN_CSF_062617_47.mzML",
        "TN_CSF_062617_49.mzML",
        "TN_CSF_062617_55.mzML",
        "TN_CSF_062617_56.mzML",
        "TN_CSF_062617_58.mzML",
        "TN_CSF_062617_60.mzML",
        "TN_CSF_062617_61.mzML",
        "TN_CSF_062617_62.mzML",
        "TN_CSF_062617_63.mzML",
    ]

    experiment_setup = {}
    i = 0
    mzml_files = []
    mzml_dir = pathlib.Path(mzml_folder)
    for file in healthy_patient_files:
        i += 1
        mzml_path = pathlib.Path(file)
        experiment_setup[str(i)] = {
            "FileName": mzml_path.stem,
            "Condition": "H",
            "Biorep": i,
            "Fraction": 1,
            "Techrep": 1,
        }
        mzml_files.append(str(mzml_dir / mzml_path))
    j = 0
    for file in als_patient_files:
        i += 1
        j += 1
        mzml_path = pathlib.Path(file)
        experiment_setup[str(i)] = {
            "FileName": mzml_path.stem,
            "Condition": "A",
            "Biorep": j,
            "Fraction": 1,
            "Techrep": 1,
        }
        mzml_files.append(str(mzml_dir / mzml_path))
    uc.params["experiment_setup"] = experiment_setup
    uc.params["quantification_evidences"] = merged_result
    quantified_peaks = uc.quantify(
        input_file=mzml_files,
        multi=True,
        engine='flash_lfq_1_1_1'
    )


if __name__ == "__main__":
    main(
        sys.argv[1],
        sys.argv[2],
    )
    pass
