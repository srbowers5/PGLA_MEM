#ifndef READ_FILE_H
#define READ_FILE_H 

class ReadFile {


    public:
        typedef struct parsed_struct_ {
            char * line;
            std::vector <char *> textFields;
            std::vector <int> intFields;
            std::vector <float> floatFields;
        } PARSED_STRUCT, * PARSED_STRUCT_PTR;

        ReadFile(char * inFile) {

            std::string line;
            int retVal;
            int i;
            int cnt;
            char * linePtr;

            std::ifstream inFd(inFile);
            if (inFd.is_open() == false) {
                printf("Failed to open %s\n", inFile);
                exit(99);
            }
            retVal = 1;
            cnt = 0;
            while (std::getline(inFd, line)) {
                linePtr = (char *) malloc(line.length() + 1);
                memcpy(linePtr, line.c_str(), line.length());
                *(linePtr+line.length()) = 0;
                lineList.push_back(linePtr);
                cnt++;
//                if (cnt >= 10000) {
//                    break;
//                }
            }
            inFd.close();
            printf("READ FILE\n");
        }


        ~ReadFile() {
             std::vector<ReadFile::PARSED_STRUCT>::iterator it2;
             std::vector<char *>::iterator it3;
             PARSED_STRUCT_PTR  parsed_line_ptr;
             char * line;
             printf("REMOVE READ FILE\n");
             while (!lineList.empty()) {
                 free(lineList.back());
                 lineList.pop_back();
             }
             while (!compSpaceList.empty()) {
                 free(compSpaceList.back());
                 compSpaceList.pop_back();
             }
             for (it2 = std::begin(parsedLine); it2 != std::end(parsedLine); ++it2) {
                 line = it2->line;
                 free(line);
                 it2->line = 0;
             }
            printf("Rm ReadFile \n");
            fflush(stdout);
        }
        void comp_space() {
            std::vector<char *>::iterator it;
            int off, off1;
            int lastSpace;
            char * currLine;
            char * linePtr;
            char nextChar;
            size_t lineLen;
            for (it = std::begin(lineList); it != std::end(lineList); ++it) {
                currLine = *it;
                linePtr = (char *) malloc(strlen(currLine) + 1);
                for (off=0; off<strlen(currLine); off++) {
                    nextChar = currLine[off];
                    if ((nextChar == '\t') || (nextChar == '\r') || (nextChar == '\n')) {
                        nextChar = ' ';
                    }
                    linePtr[off] = nextChar;
                }
                linePtr[off] = 0;
                lineLen = strlen(linePtr) + 1;
                lastSpace = 1;
                while (lineLen != strlen(linePtr)) {
                    lineLen = strlen(linePtr);
                    off1 = 0;
                    for (off=0; off<strlen(linePtr); off++) {
                        if (linePtr[off] == ' ') {
                            if (lastSpace == 1) {
                                continue;
                            } else {
                                linePtr[off1] = linePtr[off];
                                lastSpace = 1;
                                off1++;
                            }
                        } else {
                            linePtr[off1] = linePtr[off];
                            lastSpace = 0;
                            off1++;
                        }
                    }
                    if (lastSpace == 1) {
                        off1--;
                    }
                    linePtr[off1] = 0;
                }
                compSpaceList.push_back(linePtr);
            }
            printf("Comp space \n");
            fflush(stdout);
        }

        void split_lines() {
            PARSED_STRUCT parsed_str;
            std::vector<char *>::iterator it;
            PARSED_STRUCT_PTR parsed_struct_ptr;
            int off;
            int prev_was_space;
            char * linePtr;
            char * currLine;
            char * nextOff;
            int cnt;
            cnt = 0;
            for (it = std::begin(compSpaceList); it != std::end(compSpaceList); ++it) {
                parsed_str.textFields.clear();
                parsed_struct_ptr = &(parsed_str);
                currLine = *it;
                linePtr = (char *) malloc(strlen(currLine) + 1);
                if (linePtr == 0) {
                    printf("NO memory %d\n", cnt);
                    break;
                }
                memcpy(linePtr, currLine, (strlen(currLine) + 1));
                parsed_struct_ptr->line = linePtr;
                prev_was_space = 1;
                for (off = 0; off < strlen(currLine); off++) {
                     if (prev_was_space == 1) {
                         nextOff = &(linePtr[off]);
                         parsed_struct_ptr->textFields.push_back(nextOff);
                     }
                     if (linePtr[off] == ' ') {
                         prev_was_space = 1;
                         linePtr[off] = 0;
                     } else {
                         prev_was_space = 0;
                     }
                }
                parsedLine.push_back(parsed_str);
                cnt++;
                printf("Split lines %d\n", cnt);
                fflush(stdout);
            }
            printf("Split lines \n");
            fflush(stdout);
        }

        void get_float_fields() {
            std::vector<ReadFile::PARSED_STRUCT>::iterator it2;
            std::vector<char *>::iterator it3;
            float floatVal;
            for (it2 = std::begin(parsedLine); it2 != std::end(parsedLine); ++it2) {
                for (it3 = std::begin(it2->textFields); it3 != std::end(it2->textFields); it3++) {
                    floatVal = atof(*it3);
                    it2->floatFields.push_back(floatVal);
                }
            }
            printf("Got float \n");
            fflush(stdout);
        }
        std::vector <char *> lineList;
        std::vector <char *> compSpaceList;
        std::vector <PARSED_STRUCT> parsedLine;

    private:

};
#endif /* READ_FILE_H */
