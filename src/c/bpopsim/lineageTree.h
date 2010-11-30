#ifndef lineageTree_h
#define lineageTree_h

class lineageTree
{
   private:
     std::vector<long double> m_tree;
     std::vector<int> m_pointers;
   
   public:
     lineageTree()
     {
     }
     
   virtual ~lineageTree() {};

   void AddNode(long double in_mutation){ m_tree.push_back(in_mutation);}
   void SetPointer(int in_previous){m_pointers.push_back(in_previous);}
   void ClearRuns(){m_tree.clear(); m_pointers.clear();}
   void PrintTree(const std::string& output_file_name)
   {
      //Print everything out
      std::ofstream output_file;
      output_file.open(output_file_name.c_str(),std::ios_base::app);
      output_file << "tree ";
      int current_position;
      for (int on_branch=GetSizeOfTree()-1; on_branch > 0; on_branch--) 
      {
         current_position = on_branch;
         while (m_tree[current_position]!=0)
         {
            output_file << m_tree[current_position] << " ";
            current_position = m_pointers[current_position];            
         }
         output_file << std::endl;
      }
      output_file.close();
   }

   std::vector<long double> GetRuns() { return m_tree;}
   std::vector<int> GetPointers() { return m_pointers;}
   int GetSizeOfTree() {return m_tree.size();}
}; 
#endif
