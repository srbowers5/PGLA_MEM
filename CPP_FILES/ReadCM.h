
class ReadCM {

#define NUM_AAS      42
#define NUM_LIPIDS   98

    typedef struct step_cont_ {
        char    contacts[NUM_AAS][NUM_LIPIDS];
        unsigned char num_lip_cont[NUM_LIPIDS];
    } STEP_CONT, *STEP_CONT_PTR;



    public:
        std::vector <STEP_CONT_PTR> contacts_ptr;
        int numSteps;


        ReadCM(char * fileList) {

            numSteps = 0;
            std::ifstream inFd;
            std::string line;

            inFd.open(fileList,std::ifstream::binary);
            if (inFd.is_open() == false) {
                printf("Failed to open %s\n", fileList);
                exit(99);
            }

            printf("Reading %s\n", fileList);
            while (std::getline(inFd, line)) {
//                printf("Reading CM %s\n", line.c_str());
                addData(line.c_str());
            }

        }

    private:
        void addData(const char *filename);
        void addOneLine(const char * line);


};
