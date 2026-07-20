# Default data-file locations and the fetch_data!() downloader.


const DEFAULT_STARLIB_PATH = joinpath(_PACKAGE_ROOT, "data", "starlib.dat")
const LEGACY_STARLIB_PATH = joinpath(_PACKAGE_ROOT, "starlib.dat")
const DEFAULT_REACLIB_PATH = joinpath(_PACKAGE_ROOT, "data", "reaclib_v1.0.dat")
const DEFAULT_REACLIB_IL01_PATH = joinpath(_PACKAGE_ROOT, "data", "reaclib_il01.dat")
const DEFAULT_REACLIB_NACR_PATH = joinpath(_PACKAGE_ROOT, "data", "reaclib_nacr.dat")
const DEFAULT_ILIADIS2001_PATH = joinpath(_PACKAGE_ROOT, "data", "iliadis2001_rates.dat")
const DEFAULT_NACRE_PATH = joinpath(_PACKAGE_ROOT, "data", "nacre_rates.dat")
const DEFAULT_AME_PATH = joinpath(_PACKAGE_ROOT, "data", "ame2020_mass.txt")
const DEFAULT_WINVNE_PATH = joinpath(_PACKAGE_ROOT, "data", "winvne_v2.0.dat")
const STARLIB_ROWS_PER_REACTION = 60

function _default_starlib_path()
    isfile(DEFAULT_STARLIB_PATH) && return DEFAULT_STARLIB_PATH
    isfile(LEGACY_STARLIB_PATH) && return LEGACY_STARLIB_PATH
    return DEFAULT_STARLIB_PATH
end

function _default_reaclib_path()
    isfile(DEFAULT_REACLIB_PATH) || error(
        "REACLIB library not found at $(DEFAULT_REACLIB_PATH); run ReacNetJl.fetch_data!() or data/download_rates.sh to fetch it",
    )
    return DEFAULT_REACLIB_PATH
end

"""
    fetch_data!(; force=false, include_starlib=true)

Download the large rate-library data files into `data/`: the ReaclibV1.0 and
current-default JINA REACLIB snapshots, the winvne nuclear-data file, the
AME2020 atomic mass evaluation, and (unless `include_starlib=false`) the
STARLIB v6.10 library. Files already present are kept unless `force=true`.

Small provenance-critical files (the JINA `il01`/`nacr` label sets and the
tabulated paper rates) are tracked in the Git repository and never fetched.
This mirrors how pynucastro vendors its data, adapted to keep the repository
small: sources and versions are documented in `data/README.md`, and
`data/download_rates.sh` performs the same downloads from the shell.
"""
function fetch_data!(; force::Bool=false, include_starlib::Bool=true)
    data_dir = joinpath(dirname(@__DIR__), "data")
    mkpath(data_dir)

    downloads = [
        (DEFAULT_REACLIB_PATH,
         "https://reaclib.jinaweb.org/difout.php?action=cfreaclib2&library=ReaclibV1.0&rateall=1&cached=&no910=0",
         "JINA REACLIB snapshot ReaclibV1.0 (~17 MB)"),
        (joinpath(data_dir, "reaclib_default.dat"),
         "https://reaclib.jinaweb.org/difout.php?action=cfreaclib2&library=default&rateall=1&cached=&no910=0",
         "current JINA REACLIB default library (~19 MB)"),
        (DEFAULT_WINVNE_PATH,
         "https://reaclib.jinaweb.org/associated_files/Recommended/winvne_v2.0.dat",
         "winvne nuclear data (~3 MB)"),
        (DEFAULT_AME_PATH,
         "https://www-nds.iaea.org/amdc/ame2020/mass.mas20.txt",
         "AME2020 atomic mass evaluation (~0.5 MB)"),
    ]

    for (path, url, description) in downloads
        if isfile(path) && !force
            continue
        end
        @info "downloading $description"
        Downloads.download(url, path)
    end

    if include_starlib && (!isfile(DEFAULT_STARLIB_PATH) || force)
        zip_path = joinpath(data_dir, "starlib_v610.dat.zip")
        @info "downloading STARLIB v6.10 (~15 MB zip, unpacks to ~170 MB)"
        Downloads.download(
            "https://raw.githubusercontent.com/Starlib/Rate-Library/master/data/starlib_v610_120222.dat.zip",
            zip_path,
        )
        try
            run(pipeline(`unzip -o $zip_path -d $data_dir`; stdout=devnull))
        catch
            error("could not run `unzip` on $zip_path; extract it manually to $(DEFAULT_STARLIB_PATH)")
        end
        rm(zip_path; force=true)
        if !isfile(DEFAULT_STARLIB_PATH)
            extracted = filter(name -> startswith(name, "starlib") && endswith(name, ".dat"), readdir(data_dir))
            isempty(extracted) && error("STARLIB extraction did not produce a .dat file in $data_dir")
            mv(joinpath(data_dir, first(extracted)), DEFAULT_STARLIB_PATH; force=true)
        end
    end

    return nothing
end
