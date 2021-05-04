#!/usr/bin/env xcrun swift

// ----------------------------------------------------------------------
// takes a pattern and runs runTreeReader over all matching files in directory provided
//
//
// Usage:
// ../swift/runTreeReader.swift [-d ../_build] [-l] -p pattern
//
// History: 2021/05/04 First shot
//
// ----------------------------------------------------------------------

import Foundation

func getFileList(dir: String, pattern: String) -> [String] {
    var files : [String] = []

    let fileManager = FileManager.default
    do {
        let lfileList = try fileManager.contentsOfDirectory(atPath: dir)
        for file  in lfileList {
            if file.contains(pattern) && file.contains(".root") {
                files.append(file)
            }

        }
    } catch {
        print("Error while getting files: \(error.localizedDescription)")
    }

    return files
}


// ----------------------------------------------------------------------
func main() {
    // -- check that both arguments are given
    guard CommandLine.arguments.count > 2 else {
        print("usage:   \(CommandLine.arguments[0]) [-d ../_build] [-l] -p aerogel")
        print("example: ../swift/runTreeReader.swift -d ../_build/ -p cFoil -l")
        return
    }

    var logfile = false
    var pattern = "pattern"
    var dirname = "../_build"

    // -- command line parsing
    while case let option = getopt(CommandLine.argc, CommandLine.unsafeArgv, "d:lp:"), option != -1 {
        switch UnicodeScalar(CUnsignedChar(option)) {
        case "d":
            dirname = String(cString: optarg)
        case "l":
            logfile = true
        case "p":
            pattern = String(cString: optarg)
        default:
            return
        }
    }

    let files = getFileList(dir: dirname, pattern: pattern)

    for ifile in files {
        if logfile {
            let logfile = ifile.replacingOccurrences(of: ".root", with: ".log")
            print("bin/runTreeReader01 -f \(dirname)/\(ifile) >&! \(logfile)")
        } else {
            print("bin/runTreeReader01 -f \(dirname)/\(ifile) ")
        }

    }

}

// ----------------------------------------------------------------------
main()
exit(EXIT_SUCCESS)
// ----------------------------------------------------------------------
