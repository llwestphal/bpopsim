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
        m_tree.push_back(0);
	m_pointers.push_back(0);
     }
     
   virtual ~lineageTree() {};

   void AddNode(long double in_mutation){ m_tree.push_back(in_mutation);}
   void SetPointer(int in_previous){m_pointers.push_back(in_previous);}
   void ClearRuns(){m_tree.clear(); m_pointers.clear();}

   std::vector<long double> GetRuns() { return m_tree;}
   std::vector<int> GetPointers() { return m_pointers;}
}; 
#endif
