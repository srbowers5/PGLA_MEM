 
class ParseLine {

    public:
        ParseLine(char * fileName) {
            inFd = fopen(fileName, "r");
            if (inFd == 0) {
                printf("Cannot open file (%s)\n", fileName);
                exit(90);
            }
            return;
        }

        ~ParseLine() {
            fclose(inFd);
        }
        
//
//  Read one line and return a vector of field pointers
//  It will remove leading and trailing whitespace.

        int get_line_fields(std::vector <char *> * fields_ptr, char ** linePtr, size_t * lenPtr) {
            int nRead;
            int offset;
            char nextChar;
            bool found;
            char * line;

            fields_ptr->clear();
            if ((nRead = getline(linePtr, lenPtr, inFd)) <= 0) {
                printf("DONE GET_LINE_FIELDS\n");
                return -1;
            }
            offset = 0;
            found = false;
            line = *linePtr;
            while((offset < nRead) && (line[offset] != 0)) {
                nextChar = line[offset];
                if ((nextChar != ' ') && (nextChar != '\t') && (nextChar != '\r') && (nextChar != '\n')) {
                    if (found == false) {
                        fields_ptr->push_back(&(line[offset]));
//                        printf("ADD FIELDS %d off %d\n", fields_ptr->size(), offset);
                        found = true;
                    }
                } else {
                    found = false;
                    line[offset] = 0;
                }
                offset++;
            }
            return 0;
        }

    private:
        FILE * inFd;
};
