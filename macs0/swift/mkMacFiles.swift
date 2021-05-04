#!/usr/bin/env xcrun swift

// ----------------------------------------------------------------------
// takes a mac file and produces variations of it defined in a map
//
//
// Usage/examples:
// ../swift/mkMacFiles.swift -f vis.mac -n basename  -s /macs0/det/setTargetLength=5\ nm,10\ nm,15\ nm
//                          [-p /macs0/det/setTargetMaterial=Cfoil] [-p /macs0/generator/bgKinEnergy\ =26\ MeV]
//
// ../swift/mkMacFiles.swift -m 1
//
// History: 2021/04/29 First shot
//          2021/05/04 add modes to combine various runs
// ----------------------------------------------------------------------

import Foundation

// -- definition of scan and other parameter variations
let dbx = false
var scanPar = ""
var scanVals: [String] = []
var patterns: [String: String] = [:]

// ----------------------------------------------------------------------
// -- apply changes to string and write a new mac file for each
func mkMacFile(basename: String, macname: String, contents: [String]) {
    var sline = ""
    var fout = ""

    // -- loop over scan points
    for idx in scanVals {
        if dbx {print("idx: \(idx)")}
        for line in contents {
            sline = line
            for p in patterns.keys {
                if line.contains(p) {
                    sline = p.replacingOccurrences(of: " ", with: "") + " " + patterns[p]!
                    continue
                }
            }
            if line.contains(scanPar) {
                sline = scanPar + String(" ") + idx
            }
            fout +=  sline + "\n"
            //            print("\(sline)")
        }


        let cachefile =  basename + "-" + idx.replacingOccurrences(of: " ", with: "") + ".mac"
        do {
            // Write contents to file
            try fout.write(toFile: cachefile, atomically: false, encoding: String.Encoding.utf8)
        } catch{
            print("could not write cache file")
        }
        fout = ""
    }

}


// ----------------------------------------------------------------------
// -- Read in .mac file and replace all comments with nothing. Returns array of lines.
func readMacFile(_ filename: String) -> [String] {
    var contents = ""
    var result : [String] = []

    let url = URL(fileURLWithPath: filename)
    do {
        contents = try String(contentsOf: url)
    } catch {
        print("now what?")
    }

    let lines = contents.components(separatedBy: "\n")

    // -- remove all comments to trim down size
    for line in lines {
        let idx = line.firstIndex(of: "#") ?? line.endIndex
        if idx == line.startIndex {
            continue
        }
        let sline = String(line[line.startIndex..<idx])
        result.append(sline)
    }

    return result
}


// ----------------------------------------------------------------------
func main() {
    //    print("mkMacFiles   ->\(CommandLine.arguments)<-")
    // -- check that both arguments are given
    guard CommandLine.arguments.count > 2 else {
        print("usage: \(CommandLine.arguments[0]) -s /macs0/det/setTargetLength=5\\ nm,15\\ nm -p /macs0/det/setTargetMaterial=Cfoil -p /macs0/generator/bgKinEnergy=50\\ keV -n basename -f file.mac")
        return
    }

    var mode = 0
    var macname = ""
    var basename = "basename"

    // -- command line parsing
    while case let option = getopt(CommandLine.argc, CommandLine.unsafeArgv, "f:m:n:p:s:"), option != -1 {
        switch UnicodeScalar(CUnsignedChar(option)) {
        case "f":
            macname = String(cString: optarg)
        case "m":
            mode = Int(String(cString: optarg)) ?? 0
        case "n":
            basename = String(cString: optarg)
        case "p":
            let words = String(cString: optarg).components(separatedBy: "=")
            if dbx {print("-> words: \(words)")}
            if 2 == words.count {
                patterns[words[0]] = words[1]
            }
        case "s":
            let words = String(cString: optarg).components(separatedBy: "=")
            if dbx {print("-> scan: \(words)")}
            if 2 == words.count {
                let steps = words[1].components(separatedBy: ",")
                scanPar = words[0]
                scanVals += steps
            }
        default:
            return
        }
    }


    // 210504: setup for initial scan
    if 1 == mode {
        let arguments = [["-s", "/macs0/det/setTargetLength=10 um,100 um,1 mm,8 mm",
                          "-p", "/macs0/det/setTargetMaterial=Aerogel",
                          "-p", "/macs0/generator/bgKinEnergy =26 MeV",
                          "-n", "aerogel-26MeV", "-f", "vis.mac"]
                        , ["-s", "/macs0/det/setTargetLength=10 um,100 um,1 mm,8 mm",
                           "-p", "/macs0/det/setTargetMaterial=Aerogel",
                           "-p", "/macs0/generator/bgKinEnergy =5 MeV",
                           "-n", "aerogel-5MeV", "-f", "vis.mac"]
                        , ["-s", "/macs0/det/setTargetLength=10 um,100 um,1 mm,8 mm",
                           "-p", "/macs0/det/setTargetMaterial=Aerogel",
                           "-p", "/macs0/generator/bgKinEnergy =50 keV",
                           "-n", "aerogel-50keV", "-f", "vis.mac"]
                        , ["-s", "/macs0/det/setTargetLength=10 um,100 um,1 mm,8 mm",
                           "-p", "/macs0/det/setTargetMaterial=Aerogel",
                           "-p", "/macs0/generator/bgKinEnergy =500 keV",
                           "-n", "aerogel-500keV", "-f", "vis.mac"]

                        , ["-s", "/macs0/det/setTargetLength=5 nm,10 nm,15 nm",
                           "-p", "/macs0/det/setTargetMaterial=Cfoil",
                           "-p", "/macs0/generator/bgKinEnergy =5 keV",
                           "-n", "Cfoil-5keV", "-f", "vis.mac"]
                        , ["-s", "/macs0/det/setTargetLength=5 nm,10 nm,15 nm",
                           "-p", "/macs0/det/setTargetMaterial=Cfoil",
                           "-p", "/macs0/generator/bgKinEnergy =10 keV",
                           "-n", "Cfoil-10keV", "-f", "vis.mac"]
                        , ["-s", "/macs0/det/setTargetLength=5 nm,10 nm,15 nm",
                           "-p", "/macs0/det/setTargetMaterial=Cfoil",
                           "-p", "/macs0/generator/bgKinEnergy =15 keV",
                           "-n", "Cfoil-15keV", "-f", "vis.mac"]
                        , ["-s", "/macs0/det/setTargetLength=5 nm,10 nm,15 nm",
                           "-p", "/macs0/det/setTargetMaterial=Cfoil",
                           "-p", "/macs0/generator/bgKinEnergy =50 keV",
                           "-n", "Cfoil-50keV", "-f", "vis.mac"]

                        , ["-s", "/macs0/det/setTargetLength=500 nm,750 nm,1500 nm",
                           "-p", "/macs0/det/setTargetMaterial=G4_Al",
                           "-p", "/macs0/generator/bgKinEnergy =1 MeV",
                           "-n", "Al-1MeV", "-f", "vis.mac"]
                        , ["-s", "/macs0/det/setTargetLength=500 nm,750 nm,1500 nm",
                           "-p", "/macs0/det/setTargetMaterial=G4_Al",
                           "-p", "/macs0/generator/bgKinEnergy =5 MeV",
                           "-n", "Al-5MeV", "-f", "vis.mac"]
                        , ["-s", "/macs0/det/setTargetLength=500 nm,750 nm,1500 nm",
                           "-p", "/macs0/det/setTargetMaterial=G4_Al",
                           "-p", "/macs0/generator/bgKinEnergy =10 MeV",
                           "-n", "Al-10MeV", "-f", "vis.mac"]
        ]

        for args in arguments {
            let task = Process()
            task.executableURL = URL(fileURLWithPath: "/Users/ursl/fsx/geant4/geant4.10.04.p02/muamu/macs0/swift/mkMacFiles.swift")
            task.arguments = args
            do {
                try task.run()
            } catch {
                print("could not run task")
            }
        }
        return
    }

    let sfile = readMacFile(macname)
    mkMacFile(basename: basename, macname: macname, contents: sfile)

}

// ----------------------------------------------------------------------
main()
exit(EXIT_SUCCESS)
// ----------------------------------------------------------------------
