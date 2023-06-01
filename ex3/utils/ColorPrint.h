#ifndef COLORPRINT_H
#define COLORPRINT_H

#include <string>

class ColorPrint
{

private:
    // Reset
    inline static const std::string RESET = "\033[0m";  // Text Reset

    // Regular Colors
    inline static const std::string BLACK = "\033[0;30m";   // BLACK
    inline static const std::string RED = "\033[0;31m";     // RED
    inline static const std::string GREEN = "\033[0;32m";   // GREEN
    inline static const std::string YELLOW = "\033[0;33m";  // YELLOW
    inline static const std::string BLUE = "\033[0;34m";    // BLUE
    inline static const std::string PURPLE = "\033[0;35m";  // PURPLE
    inline static const std::string CYAN = "\033[0;36m";    // CYAN
    inline static const std::string WHITE = "\033[0;37m";   // WHITE
    inline static const std::string GRAY = "\033[0;90m";      // GRAY
    //const std::string MAGENTA "\033[4;35m";   // MAGENTA

    // Bold
    inline static const std::string BLACK_BOLD = "\033[1;30m";  // BLACK
    inline static const std::string RED_BOLD = "\033[1;31m";    // RED
    inline static const std::string GREEN_BOLD = "\033[1;32m";  // GREEN
    inline static const std::string YELLOW_BOLD = "\033[1;33m"; // YELLOW
    inline static const std::string BLUE_BOLD = "\033[1;34m";   // BLUE
    inline static const std::string PURPLE_BOLD = "\033[1;35m"; // PURPLE
    inline static const std::string CYAN_BOLD = "\033[1;36m";   // CYAN
    inline static const std::string WHITE_BOLD = "\033[1;37m";  // WHITE

    // Underline
    inline static const std::string BLACK_UNDERLINED = "\033[4;30m";  // BLACK
    inline static const std::string RED_UNDERLINED = "\033[4;31m";    // RED
    inline static const std::string GREEN_UNDERLINED = "\033[4;32m";  // GREEN
    inline static const std::string YELLOW_UNDERLINED = "\033[4;33m"; // YELLOW
    inline static const std::string BLUE_UNDERLINED = "\033[4;34m";   // BLUE
    inline static const std::string PURPLE_UNDERLINED = "\033[4;35m"; // PURPLE
    inline static const std::string CYAN_UNDERLINED = "\033[4;36m";   // CYAN
    inline static const std::string WHITE_UNDERLINED = "\033[4;37m";  // WHITE

    // Background
    inline static const std::string BLACK_BACKGROUND = "\033[40m";  // BLACK
    inline static const std::string RED_BACKGROUND = "\033[41m";    // RED
    inline static const std::string GREEN_BACKGROUND = "\033[42m";  // GREEN
    inline static const std::string YELLOW_BACKGROUND = "\033[43m"; // YELLOW
    inline static const std::string BLUE_BACKGROUND = "\033[44m";   // BLUE
    inline static const std::string PURPLE_BACKGROUND = "\033[45m"; // PURPLE
    inline static const std::string CYAN_BACKGROUND = "\033[46m";   // CYAN
    inline static const std::string WHITE_BACKGROUND = "\033[47m";  // WHITE

    // High Intensity
    inline static const std::string BLACK_BRIGHT = "\033[0;90m";  // BLACK
    inline static const std::string RED_BRIGHT = "\033[0;91m";    // RED
    inline static const std::string GREEN_BRIGHT = "\033[0;92m";  // GREEN
    inline static const std::string YELLOW_BRIGHT = "\033[0;93m"; // YELLOW
    inline static const std::string BLUE_BRIGHT = "\033[0;94m";   // BLUE
    inline static const std::string PURPLE_BRIGHT = "\033[0;95m"; // PURPLE
    inline static const std::string CYAN_BRIGHT = "\033[0;96m";   // CYAN
    inline static const std::string WHITE_BRIGHT = "\033[0;97m";  // WHITE

    // Bold High Intensity
    inline static const std::string BLACK_BOLD_BRIGHT = "\033[1;90m"; // BLACK
    inline static const std::string RED_BOLD_BRIGHT = "\033[1;91m";   // RED
    inline static const std::string GREEN_BOLD_BRIGHT = "\033[1;92m"; // GREEN
    inline static const std::string YELLOW_BOLD_BRIGHT = "\033[1;93m";// YELLOW
    inline static const std::string BLUE_BOLD_BRIGHT = "\033[1;94m";  // BLUE
    inline static const std::string PURPLE_BOLD_BRIGHT = "\033[1;95m";// PURPLE
    inline static const std::string CYAN_BOLD_BRIGHT = "\033[1;96m";  // CYAN
    inline static const std::string WHITE_BOLD_BRIGHT = "\033[1;97m"; // WHITE

    // High Intensity backgrounds
    inline static const std::string BLACK_BACKGROUND_BRIGHT = "\033[0;100m";// BLACK
    inline static const std::string RED_BACKGROUND_BRIGHT = "\033[0;101m";// RED
    inline static const std::string GREEN_BACKGROUND_BRIGHT = "\033[0;102m";// GREEN
    inline static const std::string YELLOW_BACKGROUND_BRIGHT = "\033[0;103m";// YELLOW
    inline static const std::string BLUE_BACKGROUND_BRIGHT = "\033[0;104m";// BLUE
    inline static const std::string PURPLE_BACKGROUND_BRIGHT = "\033[0;105m"; // PURPLE
    inline static const std::string CYAN_BACKGROUND_BRIGHT = "\033[0;106m";  // CYAN
    inline static const std::string WHITE_BACKGROUND_BRIGHT = "\033[0;107m";   // WHITE

public:
    /*
    'r' : red, 'b': blue, 'g': green, 'c': cyan, 'p': purple, 'y': yellow, 'd': darker (gray), 'n': reset colors
    */
    static std::string dye(std::string toPaint, char col) {
        switch(col) {
            case 'r':
                return RED + toPaint + RESET;
            case 'b':
                return BLUE + toPaint + RESET;
            case 'g':
                return GREEN + toPaint + RESET;
            case 'c':
                return CYAN + toPaint + RESET;
            case 'p':
                return PURPLE + toPaint + RESET;
            case 'y':
                return YELLOW + toPaint + RESET;
            case 'd':
                return GRAY + toPaint + RESET;
            default:
                return RESET + toPaint + RESET;
        }
    }
};


#endif