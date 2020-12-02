#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

def projectHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_magenta = params.monochrome_logs ? '' : "\033[0;95m";

return """  -${c_dim}-----------------------------------------------------------------${c_reset}-
${c_magenta}   ____  _   _        _____                                                          _ ${c_reset}            
${c_magenta}  / __ \\| | (_)      |  __ \\                                                        (_) ${c_reset}            
${c_magenta} | |  | | |_ _  ___  | |__) |___ _ __  _ __ ___   __ _ _ __ __ _ _ __ ___  _ __ ___  _ _ __   __ _ ${c_reset}
${c_magenta} | |  | | __| |/ __| |  _  // _ | '_ \\| '__/ _ \\ / _` | '__/ _` | '_ ` _ \\| '_ ` _ \\| | '_ \\ / _` |${c_reset}
${c_magenta} | |__| | |_| | (__  | | \\ |  __| |_) | | | (_) | (_| | | | (_| | | | | | | | | | | | | | | | (_| |${c_reset}
${c_magenta}  \\____/ \\__|_|\\___| |_|  \\_\\___| .__/|_|  \\___/ \\__, |_|  \\__,_|_| |_| |_|_| |_| |_|_|_| |_|\\__, |${c_reset}
${c_magenta}                                | |               __/ |                                       __/ |${c_reset}
${c_magenta}                                |_|              |___/                                       |___/ ${c_reset}
${c_blue}                              _             _______ _     _                                        ${c_reset}
${c_blue}                        /\\   | |           |__   __| |   (_)                                       ${c_reset}
${c_blue}                       /  \\  | | _____  __    | |  | |__  _  ___ _ __ _   _                        ${c_reset}
${c_blue}                      / /\\ \\ | |/ _ \\ \\/ /    | |  | '_ \\| |/ _ | '__| | | |                       ${c_reset}
${c_blue}                     / ____ \\| |  __/>  <     | |  | | | | |  __| |  | |_| |                       ${c_reset}
${c_blue}                    /_/    \\_|_|\\___/_/\\_\\    |_|  |_| |_|_|\\___|_|   \\__, |                       ${c_reset}
${c_blue}                                                                       __/ |                       ${c_reset}
${c_blue}                                                                      |___/                        ${c_reset}
${c_dim} alexthiery/otic-reprogramming ${c_reset}
-${c_dim}---------------------------------------------------------------${c_reset}-
    """.stripIndent()
}
