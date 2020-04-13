#include <bits/stdc++.h>
#include <cmath>
#include <chrono>
#include <thread>
#include <mutex>
#include <unistd.h>
#include <iostream>
using namespace std;


#define ROW 335
#define COL 335

// Creating a shortcut for int, int pair type
// location of the grid
typedef pair<int, int> Pair;
  
// Creating a shortcut for pair<int, pair<int, int>> type
// f_bar_1 value, location of the grid
// f_bar_2 value, location of the grid
typedef pair<double, pair<int, int>> pPair;

// A structure to hold the neccesary parameters
struct cell
{
    // Row and Column index of its parent
    // Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
    int parent_i_1, parent_j_1, parent_i_2, parent_j_2;
    // need both direction parent for back tracking as we are finding multiple paths
    // f_1 is value for f_bar_1, etc.
    double f_1, f_2, g_1, g_2, h_1, h_2;
};

ofstream output;

mutex myMutex;

char grid[ROW][COL];


Pair endpoint = make_pair(-1,-1);

// Create a closed list and initialise it to false which means
// that no cell has been included yet
// This closed list is implemented as a boolean 2D array
bool closedList[ROW][COL];

// Declare a 2D array of structure to hold the details
//of that cell
cell cellDetails[ROW][COL];


double UB = FLT_MAX; 
double optimalcost = FLT_MAX;

double gf1min;
double gf2min;
    set<pPair> openList_1;
    set<pPair> openList_2;

bool terminateF = false;
bool terminateB = false;


void printmap();
bool isValid(int row, int col);
bool isUnBlocked(int row, int col);
double calculateHValue(int row, int col, Pair dest);
void tracePath();
void expandForward(set<pPair>& openList,Pair start, Pair dest);
void expandBackward(set<pPair>& openList,Pair start, Pair dest);
void fLoop();
void bLoop();
void dibbSearch(Pair start, Pair dest);

void printmap()
{
    for(int i = 0; i < ROW; i++)
    {
      for(int j = 0; j< COL; j++)
      {
        output<< grid[i][j];
      }
      output<< endl;
    }
    output<< "\n ==================================== \n "<<endl;
}

// A Utility Function to check whether given cell (row, col)
// is a valid cell or not.
bool isValid(int row, int col)
{
    // Returns true if row number and column number
    // is in range
    return (row >= 0) && (row < ROW) &&
           (col >= 0) && (col < COL);
}
  
// A Utility Function to check whether the given cell is
// blocked or not
bool isUnBlocked(int row, int col)
{
    // Returns true if the cell is not blocked else false
    if (grid[row][col] == '#')
        return false;
    else
        return true;
}

// Manhattan Distance to calculate the 'h' heuristics.
/* input row: current position row
         col: current position col
         dest: final destination
*/
double calculateHValue(int row, int col, Pair dest)
{
    // Return using the distance formula
 // cout<<"The H value for cell "<<row<<","<<col<<" is "<<abs(((row-dest.first))
   //                       + ((col-dest.second)))<<endl;
    return ((abs(row-dest.first))
                          + (abs(col-dest.second)));
}

// A Utility Function to trace the path from the source
// to destination
void tracePath()
{
    cout<<"\nThe Path is "<<endl;

    int rowfor = endpoint.first;
    int colfor = endpoint.second;
    int rowbac = endpoint.first;
    int colbac = endpoint.second;
   // cout<< row <<" "<<col<<endl;
    stack<Pair> Pathfor;
    stack<Pair> Pathbac;
    // for(int i = 0; i < ROW; i++)
    // {
    //   for(int j = 0; j< COL; j++)
    //   {
    //     cout<<"grid:"<<i<<" ,"<<j<<"  parent: "<<cellDetails[i][j].parent_i<<" ,"<<cellDetails[i][j].parent_j<<endl;
    //   }
    //   cout<< endl;
    // }

    while (!(cellDetails[rowfor][colfor].parent_i_1 == rowfor
             && cellDetails[rowfor][colfor].parent_j_1 == colfor ))
    {
        Pathfor.push (make_pair (rowfor, colfor));
        int temp_row_1 = cellDetails[rowfor][colfor].parent_i_1;
        int temp_col_1 = cellDetails[rowfor][colfor].parent_j_1;
        rowfor = temp_row_1;
        colfor = temp_col_1;
    }

    Pathfor.push (make_pair (rowfor, colfor));

    while (!Pathfor.empty())
    {
        pair<int,int> p1 = Pathfor.top();
        Pathfor.pop();
        //printf("-> (%d,%d) ",p.first,p.second);
        grid[p1.first][p1.second]='+';
    }


    while (!(cellDetails[rowbac][colbac].parent_i_2 == rowbac
             && cellDetails[rowbac][colbac].parent_j_2 == colbac ))
    {
        Pathbac.push (make_pair (rowbac, colbac));
        int temp_row = cellDetails[rowbac][colbac].parent_i_2;
        int temp_col = cellDetails[rowbac][colbac].parent_j_2;
        rowbac = temp_row;
        colbac = temp_col;
    }
    Pathbac.push (make_pair (rowbac, colbac));
    while (!Pathbac.empty())
    {
        pair<int,int> p2 = Pathbac.top();
        Pathbac.pop();
        //printf("-> (%d,%d) ",p.first,p.second);
        grid[p2.first][p2.second]='+';
    }



    printmap();

    for(int i = 0; i < ROW; i++)
    {
      for(int j = 0; j< COL; j++)
      {
        cout<< grid[i][j];
      }
      cout<< endl;
    }
    cout<< "\n ==================================== \n "<<endl;
    return;
}


void expandForward(set<pPair>& openList,Pair start, Pair dest)
{
    pPair p1 = *openList.begin();
    double f1min = p1.first;
    int i = p1.second.first;
    int j = p1.second.second;
    
    myMutex.lock();
    openList.erase(openList.begin());
    closedList[i][j] = true;
    myMutex.unlock();

    // add i,j to closed as we are expanding all of its children
   
    /*
    Generating all the 4 successor of this cell

            N
            |
        W--Cell--E
            |
            S

    Cell-->Popped Cell (i, j)
    N -->  North       (i-1, j)
    S -->  South       (i+1, j)
    E -->  East        (i, j+1)
    W -->  West        (i, j-1)                 */

    //----------- 1st Successor (North) ------------

    // Only process this cell if this is a valid one
    if (isValid(i-1, j) == true)
    {
       
        // check if the new cell is in the open  
         if (closedList[i-1][j] == false &&
                 isUnBlocked( i-1, j) == true)
        {
            double h1new = calculateHValue(i-1,j,dest);
            double h2new = calculateHValue(i-1,j,start);
            double fbarval = 2*(cellDetails[i][j].g_1+1.0)+
                           h1new-h2new; 
            //cout << "fbarval " << fbarval <<endl;
            if(cellDetails[i-1][j].f_1>fbarval)
            {
                myMutex.lock();
                cellDetails[i-1][j].f_1=fbarval;
                cellDetails[i-1][j].g_1 = cellDetails[i][j].g_1+1.0;
                cellDetails[i-1][j].parent_i_1 = i;
                cellDetails[i-1][j].parent_j_1 = j;
                cellDetails[i-1][j].h_1 = h1new;
                cellDetails[i-1][j].h_2 = h2new;
                if(optimalcost>cellDetails[i-1][j].g_1+cellDetails[i-1][j].g_2)
                {
                    optimalcost = cellDetails[i-1][j].g_1+cellDetails[i-1][j].g_2;
                    endpoint = make_pair(i-1,j);
                    cout << optimalcost << " Optimal Cost at " << i-1 <<" ," <<j<<endl;
                }
                if(UB>(cellDetails[i-1][j].g_1+cellDetails[i-1][j].g_2))
                {
                    UB = cellDetails[i-1][j].g_1+cellDetails[i-1][j].g_2;
                    cout << "Upper Bound 1 north: " << UB << endl;
                }
                myMutex.unlock();
                
            }
          //  grid[i-1][j]='-';
          //     printmap(grid);
            openList.insert(make_pair(cellDetails[i-1][j].f_1,make_pair(i-1, j)));
        }
        
    }
    //----------- 2nd Successor (South) ------------

    if (isValid(i+1, j) == true)
    {
       
        // check if the new cell is in the open  
         if (closedList[i+1][j] == false &&
                 isUnBlocked( i+1, j) == true)
        {
            double h1new = calculateHValue(i+1,j,dest);
            double h2new = calculateHValue(i+1,j,start);
            double fbarval = 2*(cellDetails[i][j].g_1+1.0)+
                           h1new-h2new; 
            //cout << "fbarval " << fbarval <<endl;
            if(cellDetails[i+1][j].f_1>fbarval)
            {
                myMutex.lock();
                cellDetails[i+1][j].f_1=fbarval;
                cellDetails[i+1][j].g_1 = cellDetails[i][j].g_1+1.0;
                cellDetails[i+1][j].parent_i_1 = i;
                cellDetails[i+1][j].parent_j_1 = j;
                cellDetails[i+1][j].h_1 = h1new;
                cellDetails[i+1][j].h_2 = h2new;
                if(optimalcost>cellDetails[i+1][j].g_1+cellDetails[i+1][j].g_2)
                {
                    optimalcost = cellDetails[i+1][j].g_1+cellDetails[i+1][j].g_2;
                    endpoint = make_pair(i+1,j);
                    cout << optimalcost << " Optimal Cost at " << i+1 <<" ," <<j<<endl;
                }
                if(UB>(cellDetails[i+1][j].g_1+cellDetails[i+1][j].g_2))
                {
                    UB = cellDetails[i+1][j].g_1+cellDetails[i+1][j].g_2;
                    cout << "Upper Bound 1 south: " << UB << endl;
                }
                myMutex.unlock();
                
            }
          //  grid[i+1][j]='-';
          //     printmap(grid);
            openList.insert(make_pair(cellDetails[i+1][j].f_1,make_pair(i+1, j)));
        }
        
    }
 
    //----------- 3rd Successor (East) ------------
        
    if (isValid(i, j+1) == true)
    {
       
        // check if the new cell is in the open  
         if (closedList[i][j+1] == false &&
                 isUnBlocked( i, j+1) == true)
        {
            double h1new = calculateHValue(i,j+1,dest);
            double h2new = calculateHValue(i,j+1,start);
            double fbarval = 2*(cellDetails[i][j].g_1+1.0)+
                           h1new-h2new; 
            //cout << "fbarval " << fbarval <<endl;
            if(cellDetails[i][j+1].f_1>fbarval)
            {
                myMutex.lock();
                cellDetails[i][j+1].f_1=fbarval;
                cellDetails[i][j+1].g_1 = cellDetails[i][j].g_1+1.0;
                cellDetails[i][j+1].parent_i_1 = i;
                cellDetails[i][j+1].parent_j_1 = j;
                cellDetails[i][j+1].h_1 = h1new;
                cellDetails[i][j+1].h_2 = h2new;
                if(optimalcost>cellDetails[i][j+1].g_1+cellDetails[i][j+1].g_2)
                {
                    optimalcost = cellDetails[i][j+1].g_1+cellDetails[i][j+1].g_2;
                    endpoint = make_pair(i,j+1);
                    cout << optimalcost << " Optimal Cost at " << i <<" ," <<j+1<<endl;
                }                
                if(UB>(cellDetails[i][j+1].g_1+cellDetails[i][j+1].g_2))
                {
                    UB = cellDetails[i][j+1].g_1+cellDetails[i][j+1].g_2;
                    cout << "Upper Bound 1 east: " << UB << endl;
                }
                myMutex.unlock();
            }
          //  grid[i][j+1]='-';
          //     printmap(grid);
            openList.insert(make_pair(cellDetails[i][j+1].f_1,make_pair(i, j+1)));
        }
        
    }

  
    //     //----------- 4th Successor (West) ------------
    if (isValid(i, j-1) == true)
    {
       
        // check if the new cell is in the open  
         if (closedList[i][j-1] == false &&
                 isUnBlocked( i, j-1) == true)
        {
            double h1new = calculateHValue(i,j-1,dest);
            double h2new = calculateHValue(i,j-1,start);
            double fbarval = 2*(cellDetails[i][j].g_1+1.0)+
                           h1new-h2new; 
            //cout << "fbarval " << fbarval <<endl;
            if(cellDetails[i][j-1].f_1>fbarval)
            {
                myMutex.lock();
                cellDetails[i][j-1].f_1=fbarval;
                cellDetails[i][j-1].g_1 = cellDetails[i][j].g_1+1.0;
                cellDetails[i][j-1].parent_i_1 = i;
                cellDetails[i][j-1].parent_j_1 = j;
                cellDetails[i][j-1].h_1 = h1new;
                cellDetails[i][j-1].h_2 = h2new;
                if(optimalcost>cellDetails[i][j-1].g_1+cellDetails[i][j-1].g_2)
                {
                    optimalcost = cellDetails[i][j-1].g_1+cellDetails[i][j-1].g_2;
                    endpoint = make_pair(i,j-1);
                    cout << optimalcost << " Optimal Cost at " << i <<" ," <<j-1<<endl;
                }
                if(UB>(cellDetails[i][j-1].g_1+cellDetails[i][j-1].g_2))
                {
                    UB = cellDetails[i][j-1].g_1+cellDetails[i][j-1].g_2;
                    cout << "Upper Bound 1 west: " << UB << endl;
                }
                myMutex.unlock();
            }
         //   grid[i][j-1]='-';
          //     printmap(grid);
            openList.insert(make_pair(cellDetails[i][j-1].f_1,make_pair(i, j-1)));
        }
        
    }  
}

void expandBackward(set<pPair>& openList,Pair start, Pair dest)
{
    pPair p2 = *openList.begin();
    double f2min = p2.first;
    int i = p2.second.first;
    int j = p2.second.second;
    myMutex.lock();
    openList.erase(openList.begin());
    // add i,j to closed as we are expanding all of its children
    closedList[i][j] = true;
    myMutex.unlock();
    /*
    Generating all the 4 successor of this cell

            N
            |
        W--Cell--E
            |
            S

    Cell-->Popped Cell (i, j)
    N -->  North       (i-1, j)
    S -->  South       (i+1, j)
    E -->  East        (i, j+1)
    W -->  West           (i, j-1)                 */

    //----------- 1st Successor (North) ------------

    // Only process this cell if this is a valid one
    if (isValid(i-1, j) == true)
    {
       
        // check if the new cell is in the open  
         if (closedList[i-1][j] == false &&
                 isUnBlocked( i-1, j) == true)
        {
            double h1new = calculateHValue(i-1,j,dest);
            double h2new = calculateHValue(i-1,j,start);
            double fbarval = 2*(cellDetails[i][j].g_2+1.0)+
                           h2new-h1new; 

            if(cellDetails[i-1][j].f_2>fbarval)
            {
                myMutex.lock();
                cellDetails[i-1][j].f_2=fbarval;
                cellDetails[i-1][j].g_2 = cellDetails[i][j].g_2+1.0;
                cellDetails[i-1][j].parent_i_2 = i;
                cellDetails[i-1][j].parent_j_2 = j;
                cellDetails[i-1][j].h_1 = h1new;
                cellDetails[i-1][j].h_2 = h2new;
                if(optimalcost>cellDetails[i-1][j].g_1+cellDetails[i-1][j].g_2)
                {
                    optimalcost = cellDetails[i-1][j].g_1+cellDetails[i-1][j].g_2;
                    endpoint = make_pair(i-1,j);
                    cout << optimalcost << " Optimal Cost at " << i-1 <<" ," <<j<<endl;
                }
                if(UB>(cellDetails[i-1][j].g_2+cellDetails[i-1][j].g_1))
                {
                    UB = cellDetails[i-1][j].g_2+cellDetails[i-1][j].g_1;
                    cout << "Upper Bound 2 north: " << UB << endl;
                }
               myMutex.unlock();
            }
           //  grid[i-1][j]='*';
           //     printmap(grid);
            openList.insert(make_pair(cellDetails[i-1][j].f_2,make_pair(i-1, j)));
        }
        
    }
    //----------- 2nd Successor (South) ------------

    if (isValid(i+1, j) == true)
    {
       
        // check if the new cell is in the open  
         if (closedList[i+1][j] == false &&
                 isUnBlocked( i+1, j) == true)
        {
            double h1new = calculateHValue(i+1,j,dest);
            double h2new = calculateHValue(i+1,j,start);
            double fbarval = 2*(cellDetails[i][j].g_2+1.0)+
                           h2new-h1new;  

            if(cellDetails[i+1][j].f_2>fbarval)
            {
                myMutex.lock();
                cellDetails[i+1][j].f_2=fbarval;
                cellDetails[i+1][j].g_2 = cellDetails[i][j].g_2+1.0;
                cellDetails[i+1][j].parent_i_2 = i;
                cellDetails[i+1][j].parent_j_2 = j;
                cellDetails[i+1][j].h_1 = h1new;
                cellDetails[i+1][j].h_2 = h2new;
                if(optimalcost>cellDetails[i+1][j].g_1+cellDetails[i+1][j].g_2)
                {
                    optimalcost = cellDetails[i+1][j].g_1+cellDetails[i+1][j].g_2;
                    endpoint = make_pair(i+1,j);
                    cout << optimalcost << " Optimal Cost at " << i+1 <<" ," <<j<<endl;
                }                
                if(UB>(cellDetails[i+1][j].g_1+cellDetails[i+1][j].g_2))
                {
                    UB = cellDetails[i+1][j].g_1+cellDetails[i+1][j].g_2;
                    cout << "Upper Bound 2 South: " << UB << endl;
                }
                myMutex.unlock();
                
            }
           // grid[i+1][j]='*';
           //    printmap(grid);
            openList.insert(make_pair(cellDetails[i+1][j].f_2,make_pair(i+1, j)));
        }
        
    }
 
    //----------- 3rd Successor (East) ------------
        
    if (isValid(i, j+1) == true)
    {
       
        // check if the new cell is in the open  
         if (closedList[i][j+1] == false &&
                 isUnBlocked( i, j+1) == true)
        {
            double h1new = calculateHValue(i,j+1,dest);
            double h2new = calculateHValue(i,j+1,start);
            double fbarval = 2*(cellDetails[i][j].g_2+1.0)+
                           h2new-h1new;  

            if(cellDetails[i][j+1].f_2>fbarval)
            {
                myMutex.lock();
                cellDetails[i][j+1].f_2=fbarval;
                cellDetails[i][j+1].g_2 = cellDetails[i][j].g_2+1.0;
                cellDetails[i][j+1].parent_i_2 = i;
                cellDetails[i][j+1].parent_j_2 = j;
                cellDetails[i][j+1].h_1 = h1new;
                cellDetails[i][j+1].h_2 = h2new;
                if(optimalcost>cellDetails[i][j+1].g_1+cellDetails[i][j+1].g_2)
                {
                    optimalcost = cellDetails[i][j+1].g_1+cellDetails[i][j+1].g_2;
                    endpoint = make_pair(i,j+1);
                    cout << optimalcost << " Optimal Cost at " << i <<" ," <<j+1<<endl;
                }
                if(UB>(cellDetails[i][j+1].g_1+cellDetails[i][j+1].g_2))
                {
                    UB = cellDetails[i][j+1].g_1+cellDetails[i][j+1].g_2;
                    cout << "Upper Bound 2 East: " << UB << endl;
                }
                myMutex.unlock();
            }
           // grid[i][j+1]='*';
           //    printmap(grid);
            openList.insert(make_pair(cellDetails[i][j+1].f_2,make_pair(i, j+1)));
        }
        
    }

  
    //     //----------- 4th Successor (West) ------------
    if (isValid(i, j-1) == true)
    {
       
        // check if the new cell is in the open  
         if (closedList[i][j-1] == false &&
                 isUnBlocked( i, j-1) == true)
        {
            double h1new = calculateHValue(i,j-1,dest);
            double h2new = calculateHValue(i,j-1,start);
            double fbarval = 2*(cellDetails[i][j].g_2+1.0)+
                           h2new-h1new; 

            if(cellDetails[i][j-1].f_2>fbarval)
            {
                myMutex.lock();
                cellDetails[i][j-1].f_2 = fbarval;
                cellDetails[i][j-1].g_2 = cellDetails[i][j].g_2+1.0;
                cellDetails[i][j-1].parent_i_2 = i;
                cellDetails[i][j-1].parent_j_2 = j;
                cellDetails[i][j-1].h_1 = h1new;
                cellDetails[i][j-1].h_2 = h2new;
                if(optimalcost>cellDetails[i][j-1].g_1+cellDetails[i][j-1].g_2)
                {
                    optimalcost = cellDetails[i][j-1].g_1+cellDetails[i][j-1].g_2;
                    endpoint = make_pair(i,j-1);
                    cout << optimalcost << " Optimal Cost at " << i <<" ," <<j-1<<endl;
                }
                if(UB>(cellDetails[i][j-1].g_1+cellDetails[i][j-1].g_2))
                {
                    UB = cellDetails[i][j-1].g_1+cellDetails[i][j-1].g_2;
                    cout << "Upper Bound 2 West: " << UB << endl;
                }
                myMutex.unlock();
            }
           //// grid[i][j-1]='*';
            //   printmap(grid);
            openList.insert(make_pair(cellDetails[i][j-1].f_2,make_pair(i, j-1)));
        }
        
    }  
}

void dibbSearch(Pair start, Pair dest)
{

    int i, j, k ,l;
  
    for (i=0; i<ROW; i++)
    {
        for (j=0; j<COL; j++)
        {
            cellDetails[i][j].f_1 = FLT_MAX;    
            cellDetails[i][j].f_2 = FLT_MAX;
            cellDetails[i][j].g_1 = FLT_MAX;    
            cellDetails[i][j].g_2 = FLT_MAX;            
            cellDetails[i][j].parent_i_1 = -1;
            cellDetails[i][j].parent_j_1 = -1;
            cellDetails[i][j].parent_i_2 = -1;
            cellDetails[i][j].parent_j_2 = -1;
        }
    }
    
    
    // Initialising the parameters of the starting node
    i = start.first, j = start.second;
    cellDetails[i][j].f_1 = 0.0;    
    //cellDetails[i][j].f_2 = INT_MAX;
    cellDetails[i][j].g_1 = 0.0;    
    //cellDetails[i][j].g_2 = INT_MAX;
    cellDetails[i][j].h_1 = 0.0;
    //cellDetails[i][j].h_2 = calculateHValue(i,j,dest);
    cellDetails[i][j].parent_i_1 = i;
    cellDetails[i][j].parent_j_1 = j;
    
    k = dest.first, l = dest.second;
    //cellDetails[k][l].f_1 = INT_MAX;    
    cellDetails[k][l].f_2 = 0.0;
    //cellDetails[k][l].g_1 = 0.0;    
    cellDetails[k][l].g_2 = 0.0;
   // cellDetails[k][l].h_1 = calculateHValue(k,l,start);
    cellDetails[i][j].h_2 = 0.0;
    cellDetails[k][l].parent_i_2 = k;
    cellDetails[k][l].parent_j_2 = l;

  
    // Put the starting cell on the open list and set its
    // 'f' as 0
    openList_1.insert(make_pair (0.0, make_pair (i, j)));
    openList_2.insert(make_pair (0.0, make_pair (k, l)));
    memset(closedList, false, sizeof (closedList));


    pPair p1 = *openList_1.begin();
    pPair p2 = *openList_2.begin();
    gf1min = p1.first;
    gf2min = p2.first;
    int d = 1;

    // Pair lastforward;
    // Pair lastbackward;

    while(UB>((gf1min+gf2min)/2.0))
    {
        if(terminateB!=true&&terminateF!=true)
        {
            thread t1(fLoop);
            thread t2(bLoop);
            t2.join();
            t1.join();
        }
        else if(terminateB==true)
        {
            thread t1(fLoop);
            t1.join();
        }
        else if(terminateF==true)
        {
            thread t2(bLoop);
            t2.join();
        }
       // cout<< "UB: " << UB << endl; 

    }


    tracePath();
    //printmap(grid);
    return;
  
}


void fLoop()
{
    // Source is the top left-most corner
    Pair start = make_pair(1, 1);
  
    // Destination is the right-bottom most corner
    Pair dest = make_pair(ROW-2, COL-2);
  
    expandForward(openList_1,start,dest);
    pPair p1 = *openList_1.begin();
    gf1min = p1.first;
    if(gf1min == 0)
    {   
        gf1min = FLT_MAX; cout<<"Forward expansion ended."<<endl;
        myMutex.lock();
        terminateF = true;
        myMutex.unlock();
    }

}


void bLoop()
{
    // Source is the top left-most corner
    Pair start = make_pair(1, 1);
  
    // Destination is the right-bottom most corner
    Pair dest = make_pair(ROW-2, COL-2);
  
    expandBackward(openList_2,start,dest);
    pPair p2 = *openList_2.begin();
    gf2min = p2.first;
    if(gf2min==0)
    { 
        gf2min=FLT_MAX; cout<<"Backward expansion ended."<<endl;
        myMutex.lock();
        terminateB = true;
        myMutex.unlock();
    }

}

// Driver program to test above functions
int main(int argc, char** argv)
{
  ifstream file { argv[1]};
  
  //  output.open("Dibbsoutputd.out");
      output.open(argv[2]);
  if (!file.is_open()) return -1;

    /* Description of the Grid-
      --> The cell is not blocked
     #--> The cell is blocked    */
    //char grid[ROW][COL];

    for(int i = 0; i < ROW; i++)
    {
      for(int j = 0; j< COL; j++)
      {
        char read;
        file >> noskipws>>read;
        if(read!='\n'){
          grid[i][j]=read;
        //  cout<<read;
        }
        else
          file >> noskipws>>grid[i][j];
       // output<<i<<" and "<<j<<" is "<<grid[i][j]<<endl;
      
      }
      //cout<<endl;
    }
    
    
    printmap();
  
        	auto starttime = chrono::steady_clock::now();
    // Source is the top left-most corner
    Pair start = make_pair(1, 1);
  
    // Destination is the right-bottom most corner
    Pair dest = make_pair(ROW-2, COL-2);
  
    dibbSearch(start, dest);

            auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds : " 
		<< chrono::duration_cast<chrono::milliseconds>(end - starttime).count()
		<< " ms" << endl;

    output.close();
  
    return(0);
}