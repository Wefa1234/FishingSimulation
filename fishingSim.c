#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define DEBUG

/*
TO DO:
- add time logs
-extend forbiddencells
-catch exception if boat cant move
- add wave propagation -> Left border cells use sine function to get a new wave 
   value which is then allways passed to the rigth neighbor.
   Add Offset to get a diagonal wave going through map.
   Add random number for different wave hights.
   if wavehight to high then cell gets blocked -> added to forbiddenCells for boats
*/

#define SIZE 36
#define DIMX 6
#define DIMY 6
#define FRAC 0.25
#define CAPACITY 3
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

typedef enum
{
    water,
    land,
    fish1,
    fish2,
    boat1,
    boat2,
    fishfish
} gridType;

void printMap(const gridType *map)
{
    for (int i = 0; i < SIZE; i++)
    {
        switch (map[i])
        {
        case 0:
            printf("~~");
            break;
        case 1:
            printf("⚓");
            break;
        case 2:
            printf("🐟");
            break;
        case 3:
            printf("🐳");
            break;
        case 4:
            printf("🚤");
            break;
        case 5:
            printf("🦆");
            break;
        case 6:
            printf("🐟");
            break;
        }
        if (i != 0 && i % DIMX == DIMX - 1)
            printf("\n");
    }
    printf("\n\n");
}
void writeMap(char *fileName, MPI_File fh ,const gridType *map, double *timer, int step)
{
    char tile[36];
    char header[128];
    double current_t = MPI_Wtime();
    MPI_Status status;

    snprintf(header,128,"Iteration step: %03d, time: %05.2fs   |\n", step,current_t-*timer);
    MPI_File_write(fh, header,strlen(header), MPI_CHAR,&status);
    for (int i = 0; i < SIZE; i++)
    {
        switch (map[i])
        {
        case 0:
            snprintf(tile,36,"~~");
            break;
        case 1:
            snprintf(tile,36,"⚓");
            break;
        case 2:
            snprintf(tile,36,"🐟");
            break;
        case 3:
            snprintf(tile,36,"🐳");
            break;
        case 4:
            snprintf(tile,36,"🚤");
            break;
        case 5:
            snprintf(tile,36,"🦆");
            break;
        case 6:
            snprintf(tile,36,"🐟");
            break;
        }
        MPI_File_write(fh, tile, strlen(tile), MPI_CHAR,&status);

        if (i != 0 && i % DIMX == DIMX - 1){
            snprintf(tile,36,"                        |\n");
            MPI_File_write(fh, tile, strlen(tile), MPI_CHAR,&status);
        }
    }
    snprintf(header,36,"\n\n");
    MPI_File_write(fh, header,strlen(header), MPI_CHAR,&status);
}

void printArray(const gridType *map)
{
    for (int i = 0; i < SIZE; i++)
    {
        printf("%d", map[i]);
        if (i != 0 && i % DIMX == DIMX - 1)
            printf("\n");
    }
    printf("\n");
    printf("Legend: 1:land, 2:fish1, 3:fish2, 4:boat1, 5:boat2\n");
}

void printDoubleArray(const double *map)
{
    for (int i = 0; i < SIZE; i++)
    {
        printf("%.2f  ", map[i]);
        if (i != 0 && i % DIMX == DIMX - 1)
            printf("\n");
    }
    printf("\n");
}

void setUpMap(gridType *map)
{
    // randomly assign position of harbor
    int harborX = rand() % DIMX;
    int harborY = rand() % DIMY;
    map[DIMX * harborX + harborY] = land;
    // randomly assign position of fish and boat

    int count = 0;
    int posX, posY, idx;
    while (count < 4)
    {
        posX = rand() % DIMX;
        posY = rand() % DIMY;
        idx = DIMX * posX + posY;
        if (map[idx] == water)
        {
            switch (count)
            {
            case 0:
                map[idx] = fish1;
                break;
            case 1:
                map[idx] = fish2;
                break;
            case 2:
                map[idx] = boat1;
                break;
            case 3:
                map[idx] = boat2;
                break;
            default:
                break;
            }
            count++;
        }
    }
    //printArray(map);
}

int contains(const int *array, const int num)
{
    size_t n = 4; // sizeof(array) / sizeof(array[0]);
    int contain = 0;
    for (int i = 0; i < n; i++)
    {
        if (array[i] == num)
            contain = 1;
    }
    // printf("%d",contain);
    return contain;
}

void getNeighbors(const int input, int *neighbors)
{
    neighbors[0] = (input - DIMX + SIZE) % SIZE;                                      // up
    neighbors[1] = input + 1 - ((input % DIMX) / (DIMX - 1)) * DIMX;                  // right
    neighbors[2] = (input + DIMX) % SIZE;                                             // down
    neighbors[3] = (input - 1) + (((input + (DIMX - 1)) % DIMX) / (DIMX - 1)) * DIMX; // left
    // printf("%d ",neighbors[0]);
    // printf("%d ",neighbors[1]);
    // printf("%d ",neighbors[2]);
    // printf("%d ",neighbors[3]);
}

void getBlockedCells(int boat1, int boat2, int *forbiddenPos)
{
    int neighbors1[4];
    getNeighbors(boat1, neighbors1);
    int neighbors2[4];
    getNeighbors(boat2, neighbors2);
    int count = 0;
    for (int j = 0; j < 4; j++)
    {
        if (contains(neighbors2, neighbors1[j]))
        {
            forbiddenPos[count] = neighbors1[j];
            count++;
        }
    }
}

void convertIdx2Coords(const int idx1d, int *coords){
    coords[0] = idx1d / DIMX;
    coords[1] = idx1d % DIMX;
}
int convertCoords2Idx(const int coords[]){
    int idx1d = DIMX * coords[0] + coords[1];
    return idx1d;
}

int moveBoat(const int *nbrs, const int *forbiddenPos, const int harborPos)
{
    int free = 0;
    int idx;
    while (!free)
    {
        idx = rand() % 4;                                                 // pick random indx
        if (!contains(forbiddenPos, nbrs[idx]) && nbrs[idx] != harborPos) // TODO: ADD harbor avoidance
        {
            free = 1;
        }
    }
    return nbrs[idx];
}

int moveBoat2Harbor(const int *nbrs, const int *forbiddenPos, const int harborPos, const int rank,int *fishCount, double *timer, MPI_File fh, MPI_Status status, int step)
{  
    int free = 0;
    int idx;
    int count = 0;
    int dir[2];
    int harbor[2];
    int boat[2];
    char message[128];
    convertIdx2Coords(harborPos, harbor);  // get coords of harbor
    convertIdx2Coords(rank, boat);  //get coords of boat
    dir[0] = harbor[0]-boat[0];   // calculate direction to harbor
    dir[1] = harbor[1]-boat[1];
    if(dir[0]*dir[0]+dir[1]*dir[1]==1){
        (*fishCount) = 0;
        printf("setting fishcount to 0\n");

        //calculate fishing time 
        double current_t = MPI_Wtime();      
        //printf("Timer from function %1.2f.\n",*timer); //test printing  
	//printf("Current time from fuction %1.2f. \n",current_t);
        printf("Boat took %1.2f seconds to fill the nets and return to the harbor.\n",current_t-*timer);
	snprintf(message,128,"Boat arrived to harbor.        Iteration step: %03d. \n", step);
        MPI_File_write_shared(fh, message, strlen(message), MPI_CHAR, &status);
        (*timer) = MPI_Wtime(); //reset timer 
                
    }
    //printf("dir %d,%d\n", dir[0],dir[1]);
    while (!free)
    {
        if(dir[0] != 0 && count == 0){  //first try to move the boat to the row of the harbor
            if(dir[0]<0){
                idx = 0; // go up
            }else{
                idx = 1; // go down
            }
        }else if(dir[1] != 0 && (count == 1 || count == 0)){ // then try to move the boat to the column of the harbor
            if(dir[1]<0){
                idx = 2; // go left
            }else{
                idx = 3; // go right
            }
        }else{  // if both movements towards harbor are blocked pick random dir
            idx = rand() % 4;
        }
       //printf("idx = %d \n",idx);
       //printf("target cell %d\n",nbrs[idx]);
        if (!contains(forbiddenPos, nbrs[idx]) && nbrs[idx] != harborPos) // TODO: ADD harbor avoidance
        {
            free = 1;
        }else{
            count++;
        }
    }
    return nbrs[idx];
}

int moveFish(const int *nbrs, const int harborPos)
{
    int free = 0;
    int idx;
    while (!free)
    {
        idx = rand() % 4; // pick random indx
        if (nbrs[idx] != harborPos)
        {
            free = 1;
        }
    }
    return nbrs[idx];
}

int getObjPos(const gridType *map, const int type)
{
    int pos = -1;
    for (int i = 0; i < SIZE; i++)
    {
        if (map[i] == type)
        {
            pos = i;
        }
    }
    return pos;
}

void logger(MPI_Info info, char *fileName, MPI_File fh, char *text){
    MPI_Status status;
    MPI_Info_create(&info);
    //printf("Func Logger called\n");
    MPI_File_open(MPI_COMM_WORLD, fileName,MPI_MODE_CREATE | MPI_MODE_RDWR,info,&fh);
    char buf[42];
    time_t t = t;
    t = time(NULL);
    snprintf(buf,42,"%s: %s \n",asctime( localtime(&t)),text);
    MPI_File_write(fh, buf,strlen(buf), MPI_CHAR,&status);
    MPI_File_close(&fh);
}

double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

double generateWave(const int timestep,const int rank){
    double height;
    double a = 1;
    
    double PI2 = 2*M_PI;
    double k = PI2/DIMX;
    double omega = PI2/DIMX;
    //double k = 2*pi/lambda; //lambda is wavelength in m
    //double omega = 2*pi/T // T is period in s
    height = a*sin(k*timestep+rank)+a+randfrom(0,0.2);
    return height;
}

void updateState(int *inbuf, int *state, int *gotFish){
    if ((*state) != land)
        {
            (*gotFish) = 0;
            (*state) = water;
            int i = 0;
            while((*state) == water && i<4){
                switch(inbuf[i]){
                    case boat1: (*state) = boat1; break;
                    case boat2: (*state) = boat2; break;
                    case fish1: (*state) = fish1; break;
                    case fish2: (*state) = fish2; break;
                    default: break;
                }
                i++;
            }
            // In case multiple objects recv
            for(i=i;i<4;i++){
                if(inbuf[i]!= water){
                   //printf("Multiple states recv \n");
                   //printf("rank= %d with coords= (%d,%d) has received the following states from its neighbors (u,d,l,r)= %d %d %d %d \n", rank,coords[0],coords[1], inbuf[0], inbuf[1], inbuf[2], inbuf[3] );
                    if (((*state) == fish1 || (*state)==fish2) && (inbuf[i] == fish1 || inbuf[i]==fish2) ){
                        (*state) = fishfish;
                        //printf("Went to case recv second fish\n");
                    }
                    if(((*state) == fish1 || (*state) == fish2 || (*state) == fishfish) && inbuf[i]>fish2){
                        (*gotFish) = (*state);
                        (*state) = inbuf[i];
                        //printf("Went to case recv boat after fish\n");
                    }
                    if(((*state) == boat1 || (*state) == boat2) && inbuf[i]<boat1){
                        if((*gotFish) != 0){
                            (*gotFish) = fishfish;
                            //printf("Went to case recv fish after boat and had fish already\n");
                        }else{
                            (*gotFish) = inbuf[i];
                            //printf("Went to case recv fish after boat\n");
                        }
                    }
                }
                
            }
        }
}

int main(int argc, char **argv)
{
    int numtasks, rank, source, dest, i, tag = 1;
    int inbuf[4] = {MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL};
    int outbuf[4] = {MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL};
    int boatPos[2], boatPosNew[2];
    int fish1Dest;
    int fish2Dest;
    int forbiddenPos[4] = {-1, -1, -1, -1};
    int nbrs[4];
    int dims[2] = {DIMX, DIMY}, periods[2] = {1, 1}, reorder = 1;
    int coords[2];
    MPI_Comm cartcomm;

    gridType map[SIZE] = {water}; // create map variable that stores all states and can be printed
    int idx1d;
    //int free;
    int recvBuff[SIZE];
    int state;
    int helper;
    int helper2;
    int gotFish = 0;
    int fishCount[2] = {0};
    int info[2];
    const int netCapacity = CAPACITY;
    int* pointer;
    int* statePointer;
    int* gotFishPointer;

    // vars for wave propagation
    double waveHeight, waveHeightNew;
    double waveMap[SIZE] = {0};
    
    // timer variables
    double timer_boat1, timer_boat2, timer_global;

    MPI_Request reqs[8];
    //MPI_Request reqsBoats;
    //MPI_Request reqsHarbor[2];
    MPI_Status stats[8];
    MPI_Status statusBoats;
    MPI_Status statusBoats2;
    MPI_Status statusHarbor[2];

    
    // starting with MPI program
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);

    MPI_Comm_rank(cartcomm, &rank);

    // MPI_Cart_rank()
    MPI_Cart_coords(cartcomm, rank, 2, coords);

    srand(time(NULL)+rank);
    //---------insert code here-------------
    // root process creates map
    if (rank == 0)
    {
        setUpMap(map);
    }
    // Root Process broadcasts map
    MPI_Bcast(map, SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    int harborPos;
    harborPos = getObjPos(map, 1);

    // int harborPos = 10;
    //  get own idx
    idx1d = DIMX * coords[0] + coords[1];
    state = map[idx1d];
    if (state == land)
    {
        boatPos[0] = getObjPos(map, 4);
        boatPos[1] = getObjPos(map, 5);
    }
    // get their neighboring cells
    MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
    MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);

    //start fishing timers and global timer
    timer_boat1 = MPI_Wtime();
    timer_boat2 = MPI_Wtime();
    timer_global = MPI_Wtime();
    //initialize MPI i/o
    MPI_File fhMap;
    MPI_Info fileInfoMap;
    MPI_Status statusMap;
    MPI_File fhParallel;
    MPI_Info fileInfoParallel;
    MPI_Status statusParallel;
    char *fileNameMap = "MapLog";
    char *fileNameParallel = "ParallelLog";
    char eventMessage[128];
    MPI_Info_create(&fileInfoMap);
    MPI_File_open(MPI_COMM_WORLD, fileNameMap, MPI_MODE_CREATE | MPI_MODE_RDWR, fileInfoMap, &fhMap);

    MPI_Info_create(&fileInfoParallel);
    MPI_File_open(MPI_COMM_WORLD, fileNameParallel, MPI_MODE_CREATE | MPI_MODE_RDWR, fileInfoParallel, &fhParallel);


    for (int timeStep = 0; timeStep < 100; timeStep++)
    {
        outbuf[0] = 0;
        outbuf[1] = 0;
        outbuf[2] = 0;
        outbuf[3] = 0;
        if(state != land){
            boatPosNew[0] = -1;
            boatPosNew[1] = -1;
        }
        if(state == fishfish){
            //printf("FISH FISH!\n");
        }
        if (state == land)
        {
            // Harbor needs to send fishCount and pos of other boat
            info[0] = boatPos[1];
            info[1] = fishCount[0];
            MPI_Send(info, 2, MPI_INT, boatPos[0], tag, MPI_COMM_WORLD);
            info[0] = boatPos[0];
            info[1] = fishCount[1];
            MPI_Send(info, 2, MPI_INT, boatPos[1], tag, MPI_COMM_WORLD);
            // calculate forbidden pos
            forbiddenPos[0] = -1;
            forbiddenPos[1] = -1;
            forbiddenPos[2] = -1;
            forbiddenPos[3] = -1;
            //printf("Calculating blocked cells\n");
            getBlockedCells(boatPos[0], boatPos[1], forbiddenPos);
            //printf("Blocked Cells are %d,%d,%d,%d \n", forbiddenPos[0], forbiddenPos[1], forbiddenPos[2], forbiddenPos[3]);

            //---Process Waves----
            for(int j=0; j<SIZE;j++){
                if(waveMap[j]>2.15){
                    printf("Wave >2.15 encounterd blocking Cell\n");
                    int found = 0;
                    for(int k=0;k<3;k++){ //only add stormcells to forbidden list until 3 entries reached -> leave a space to go for ship TODO:
                        if(!found && forbiddenPos[k]==-1){
                            forbiddenPos[k] = j;
                            printf("added stormcells\n");
                            found = 1;
			    //write event to log file
			    snprintf(eventMessage,128,"Storm detected!                Iteration step: %03d.\n",timeStep);
                            MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
                        }else if(k==2 && !found){
                            printf("didnt add storm cell as to not block boat\n");
                        }
                    }
                }
            }

            // sende aktuell verbotene felder oder -1
            // printf("trying to send msg to %d and %d (boats)", boatPos[0],boatPos[1]);
            MPI_Send(forbiddenPos, 4, MPI_INT, boatPos[0], tag, MPI_COMM_WORLD); //, &reqsHarbor[0]);
            MPI_Send(forbiddenPos, 4, MPI_INT, boatPos[1], tag, MPI_COMM_WORLD); //, &reqsHarbor[1]);

            // recv new target pos boats
            // printf("trying to rcv msgs from %d and %d (boats)", boatPos[0],boatPos[1]);
            MPI_Recv(&helper, 1, MPI_INT, boatPos[0], tag, MPI_COMM_WORLD, &statusHarbor[0]); //, &reqsHarbor[0]);
            boatPosNew[0] = helper;
            MPI_Recv(&helper, 1, MPI_INT, boatPos[1], tag, MPI_COMM_WORLD, &statusHarbor[1]); //, &reqsHarbor[1]);
            boatPosNew[1] = helper;

            
            //recv new fishcounts
            MPI_Recv(&helper, 1, MPI_INT, boatPos[0], tag, MPI_COMM_WORLD, &statusHarbor[0]);
            fishCount[0] = helper;
            MPI_Recv(&helper, 1, MPI_INT, boatPos[1], tag, MPI_COMM_WORLD, &statusHarbor[1]);
            fishCount[1] = helper;

            boatPos[0]=boatPosNew[0];
            boatPos[1]=boatPosNew[1];

        }
        if (state == boat1)
        {
            // recv msg
            MPI_Recv(info, 2, MPI_INT, harborPos, tag, MPI_COMM_WORLD, &statusBoats);
            // printf("trying to recv msg from %d (harbor)", harborPos);
            MPI_Recv(forbiddenPos, 4, MPI_INT, harborPos, tag, MPI_COMM_WORLD, &statusBoats); //, &reqsBoats[0]);
            // calculate new target pos
            if(info[1]<netCapacity){
                boatPosNew[0] = moveBoat(nbrs, forbiddenPos, harborPos);
            }else{
                printf("Moving Boat 1 to Harbor\n");
		//write event to file
		snprintf(eventMessage,128,"Moving boat 1 to harbor.       Iteration step: %03d.\n",timeStep);
                MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
		
                pointer = &info[1];
                boatPosNew[0] = moveBoat2Harbor(nbrs, forbiddenPos, harborPos, rank, pointer, &timer_boat1, fhParallel, statusParallel ,timeStep);
            }

            for (i = 0; i < 4; i++)
            {
                if (nbrs[i] == boatPosNew[0])
                {
                    outbuf[i] = boat1; 
                }
            }
            // send new target pos to harbor
            helper = boatPosNew[0];
            // printf("trying to recv msg from %d (harbor)", harborPos);
            MPI_Send(&helper, 1, MPI_INT, harborPos, tag, MPI_COMM_WORLD);

            
            //printf("Recieved following info:%d,%d\n",info[0],info[1]);
            //check if caught fish
            if(info[1]<netCapacity){
                switch(gotFish){
                    case fish1:
                        info[1]++;
                        snprintf(eventMessage,128,"Boat 1 here, we caught %d fish. Iteration step: %03d.\n",info[1], timeStep);
                        MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
                        break;

                    case fish2:
                        info[1]++;
                        snprintf(eventMessage,128,"Boat 1 here, we caught %d fish. Iteration step: %03d.\n",info[1], timeStep);
                        MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
                        break;

                    case fishfish:
                         info[1] += 2;
                         snprintf(eventMessage,128,"Boat 1 here, we caught %d fish. Iteration step: %03d.\n",info[1], timeStep);
                         MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
                         break;
                    default: break;
                }
            }
            //send new fishcount back
            helper = info[1];
            MPI_Send(&helper, 1, MPI_INT, harborPos, tag, MPI_COMM_WORLD);

            // Chatter with second boat
            MPI_Sendrecv(&helper,1,MPI_INT,info[0],1,&helper2,1,MPI_INT,info[0],1,MPI_COMM_WORLD,&statusBoats);
            printf("Over Radio: Boat 2 here, we caught %d fish\n",helper2);
        }
        if (state == boat2)
        {
            MPI_Recv(info, 2, MPI_INT, harborPos, tag, MPI_COMM_WORLD, &statusBoats);
            // recv msg
            MPI_Recv(forbiddenPos, 4, MPI_INT, harborPos, tag, MPI_COMM_WORLD, &statusBoats); //, &reqsBoats[0]); //cant use Irec/Isend since we need forbiddenPos before we continue
            // calculate new target pos

            if(info[1]<netCapacity){
                boatPosNew[1] = moveBoat(nbrs, forbiddenPos, harborPos);
            }else{
                printf("Moving Boat 2 to Harbor\n");
		// write event to file
		snprintf(eventMessage,128,"Moving boat 2 to harbor.       Iteration step: %03d.\n",timeStep);
                MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
                
                pointer = &info[1];
                boatPosNew[1] = moveBoat2Harbor(nbrs, forbiddenPos, harborPos, rank, pointer, &timer_boat2, fhParallel, statusParallel,timeStep);
                
            }            

            for (i = 0; i < 4; i++)
            {
                if (nbrs[i] == boatPosNew[1])
                {
                    outbuf[i] = boat2; // boatPosNew[1];
                }
            }
            // send new target pos to harbor
            helper = boatPosNew[1];
            MPI_Send(&helper, 1, MPI_INT, harborPos, tag, MPI_COMM_WORLD);
            
            
            //printf("Recieved following info:%d,%d\n",info[0],info[1]);
            if(info[1]<netCapacity){
                switch(gotFish){
                    case fish1:
                        info[1]++;
                        snprintf(eventMessage,128,"Boat 2 here, we caught %d fish. Iteration step: %03d.\n",info[1], timeStep);
                        MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
                        break;

                    case fish2:
                        info[1]++;
                        snprintf(eventMessage,128,"Boat 2 here, we caught %d fish. Iteration step: %03d.\n",info[1], timeStep);
                        MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
                        break;

                    case fishfish:
                         info[1] += 2;
                         snprintf(eventMessage,128,"Boat 2 here, we caught %d fish. Iteration step: %03d.\n",info[1], timeStep);
                         MPI_File_write_shared(fhParallel, eventMessage, strlen(eventMessage), MPI_CHAR, &statusParallel);
                         break;
                    default: break;
                }
            }
            //send new fishcount back
            helper = info[1];
            MPI_Send(&helper, 1, MPI_INT, harborPos, tag, MPI_COMM_WORLD);

            MPI_Sendrecv(&helper,1,MPI_INT,info[0],1,&helper2,1,MPI_INT,info[0],1,MPI_COMM_WORLD,&statusBoats);
            printf("Over Radio: Boat 1 here, we caught %d fish\n",helper2);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (state == fish1 || gotFish == fish1)
        {
            do{
                fish1Dest = moveFish(nbrs, harborPos);
            }while(boatPosNew[0] == fish1Dest || boatPosNew[1] == fish1Dest);
            
            for (i = 0; i < 4; i++)
            {
                if (nbrs[i] == fish1Dest)
                {
                    outbuf[i] = fish1; // fish1Dest;
                }
            }
        }
        if (state == fish2 || gotFish == fish2)
        {
            do{
                fish2Dest = moveFish(nbrs, harborPos);
            }while(boatPosNew[0] == fish2Dest || boatPosNew[1] == fish2Dest);
            
            for (i = 0; i < 4; i++)
            {
                if (nbrs[i] == fish2Dest)
                {
                    outbuf[i] = fish2; // fish2Dest;
                }
            }
        }
        if (state == fishfish || gotFish == fishfish)
        {
            do{
                fish1Dest = moveFish(nbrs, harborPos);
            }while(boatPosNew[0] == fish1Dest || boatPosNew[1] == fish1Dest);
            //printf("I am fish fish\n");
            
            do{
                fish2Dest = moveFish(nbrs, harborPos);
            }while (fish2Dest == fish1Dest || (boatPosNew[0] == fish2Dest || boatPosNew[1] == fish2Dest));
            
            for (i = 0; i < 4; i++)
            {
                if (nbrs[i] == fish1Dest)
                {
                    outbuf[i] = fish1; // fish1Dest;
                }
                if (nbrs[i] == fish2Dest)
                {
                    outbuf[i] = fish2; // fish2Dest;
                }
            }
            
        }
        // Send and recv msg from Neighbors using non-Blocking communication       
        for (i = 0; i < 4; i++)
        {
            dest = nbrs[i];
            source = nbrs[i];
            MPI_Isend(&outbuf[i], 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &reqs[i]);
            MPI_Irecv(&inbuf[i], 4, MPI_INT, source, tag, MPI_COMM_WORLD, &reqs[i + 4]);
        }
        MPI_Waitall(8, reqs, stats);
        // printf("rank= %d with coords= (%d,%d) has sent the following states to its neighbors (u,d,l,r)= %d %d %d %d \n", rank,coords[0],coords[1], outbuf[0], outbuf[1], outbuf[2], outbuf[3] );
        // printf("rank= %d with coords= (%d,%d) has received the following states from its neighbors (u,d,l,r)= %d %d %d %d \n", rank,coords[0],coords[1], inbuf[0], inbuf[1], inbuf[2], inbuf[3] );
        
        //----Updating States----
        gotFishPointer = &gotFish;
        statePointer = &state;
        updateState(inbuf, statePointer, gotFishPointer);

        // ----Wave propagation-----

        if(rank%DIMX == 0){ // if left border of map generate wave that travels through map
            waveHeight = generateWave(timeStep,rank);
        }
        
        MPI_Sendrecv(&waveHeight,1,MPI_DOUBLE,nbrs[3],1,&waveHeightNew,1,MPI_DOUBLE,nbrs[2],1,MPI_COMM_WORLD,&statusBoats);
        
        waveHeight = waveHeightNew;
        
        //Lighthouse observes waves // uses gather to get wave hights -> add values of threshhold to forbidden cells;
        MPI_Gather(&waveHeight,1,MPI_DOUBLE,&waveMap,1,MPI_DOUBLE,harborPos,MPI_COMM_WORLD);
        

        //------Wait for all processes and then print the map-----
        MPI_Barrier(MPI_COMM_WORLD);
        
        MPI_Gather(&state, 1, MPI_INT, &recvBuff, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            writeMap(fileNameMap, fhMap, recvBuff, &timer_global, timeStep);
            printMap(recvBuff);
        }
        // UNCOMMENT TO SEE WAVES
        if(rank == harborPos){
            printDoubleArray(waveMap);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_File_close(&fhMap);
    MPI_File_close(&fhParallel);
    MPI_Finalize();
    return 0;
}
