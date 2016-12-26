#ifndef MEMORY_DELEGATE_H
#define MEMORY_DELEGATE_H

template<int rows, int cols, class ElemT, class MemT> class Matrix;

///////////////////////////////////////////////////////////////// Array Memory Delegate ///////////////////////////////////////////////////////////////////

template<int rows, int cols = 1, class ElemT = float> struct Array
{
    mutable ElemT m[rows * cols];

    ElemT &operator()(int row, int col) const
    {
        static ElemT dummy;

        if(row > rows || col > cols)
            return dummy;
        else
            return m[row * cols + col];
    }
};

template<int rows, int cols, class ElemT, class opElemT, class retElemT>
Matrix<rows,cols,retElemT,Array<rows,cols,retElemT> > &Add(const Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > &A, const Matrix<rows,cols,opElemT,Array<rows,cols,opElemT> > &B, Matrix<rows,cols,retElemT,Array<rows,cols,retElemT> > &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C.delegate.m[i * cols + j] = A.delegate.m[i * cols + j] + B.delegate.m[i * cols + j];

    return C;
}

template<int rows, int cols, class ElemT, class opElemT, class retElemT>
Matrix<rows,cols,retElemT,Array<rows,cols,retElemT> > &Subtract(const Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > &A, const Matrix<rows,cols,opElemT,Array<rows,cols,opElemT> > &B, Matrix<rows,cols,retElemT,Array<rows,cols,retElemT> > &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C.delegate.m[i * cols + j] = A.delegate.m[i * cols + j] - B.delegate.m[i * cols + j];

    return C;
}

template<int rows, int cols, int operandCols, class ElemT, class opElemT, class retElemT>
Matrix<rows,operandCols,retElemT,Array<rows,operandCols,retElemT> > &Multiply(const Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > &A, const Matrix<cols,operandCols,opElemT,Array<cols,operandCols,opElemT> > &B, Matrix<rows,operandCols,retElemT,Array<rows,operandCols,retElemT> > &C)
{
    int i,j,k;

    for(i = 0; i < rows; i++)
        for(j = 0; j < operandCols; j++)
        {
            if(cols > 0)
                C.delegate.m[i * operandCols + j] = A.delegate.m[i * cols] * B.delegate.m[j];

            for(k = 1; k < cols; k++)
                C.delegate.m[i * operandCols + j] += A.delegate.m[i * cols + k] * B.delegate.m[k * operandCols + j];
        }

    return C;
}

///////////////////////////////////////////////////////////////// Reference Memory Delegate ///////////////////////////////////////////////////////////////////

template<class ElemT, class MemT> struct Reference
{
    const MemT &parent;
    int rowOffset, colOffset;

    Reference<ElemT,MemT>(const MemT &obj, int rowOff, int colOff) : parent(obj), rowOffset(rowOff), colOffset(colOff) { }
    Reference<ElemT,MemT>(const Reference<ElemT,MemT> &obj) : parent(obj.parent), rowOffset(obj.rowOffset), colOffset(obj.colOffset) { }

    ElemT &operator()(int row, int col) const
    {
        return parent(row+rowOffset,col+colOffset);
    }
};

template<int rows, int cols, class ElemT = float> using ArrayRef = Reference<ElemT,Array<rows,cols,ElemT> >;
template<int rows, int cols, class ParentMemT, class ElemT = float> using RefMatrix = Matrix<rows, cols, ElemT, Reference<ElemT,ParentMemT> >;

///////////////////////////////////////////////////////////////// Identity Memory Delegate ///////////////////////////////////////////////////////////////////

template<class ElemT> struct Iden
{
    ElemT &operator()(int row, int col) const
    {
        static ElemT ret;

        if(row == col)
            return (ret = 1);
        else
            return (ret = 0);
    }
};

template<int rows, int cols, class ElemT = float> using Identity = Matrix<rows, cols, ElemT, Iden<ElemT> >;

///////////////////////////////////////////////////////////////// Zeros Memory Delegate ///////////////////////////////////////////////////////////////////

template<class ElemT> struct Zero
{
    ElemT &operator()(int row, int col) const
    {
        static ElemT ret;

        return (ret = 0);
    }
};

template<int rows, int cols = 1, class ElemT = float> using Zeros = Matrix<rows, cols, ElemT, Zero<ElemT> >;


///////////////////////////////////////////////////////////////// Sparse Memory Delegate ///////////////////////////////////////////////////////////////////

// This uses a hash table to look up row/col/val items. It uses an open addressing collision strategy so we can avoid using dynamic memory
template<int cols, int tableSize, class ElemT> struct Sparse
{
    struct HashItem
    {
        mutable int key;
        mutable ElemT val;

        HashItem() { key = -1; }

    } table[tableSize];

    ElemT &operator()(int row, int col) const
    {
        // Make a key out of the row / column
        int key = row * cols + col;

        // Calculate the hash by taking the modulo of the key with the tableSize
        int hash = key % tableSize;

        const HashItem *item;

        // Find a item with a key matching the input key
        for(int i = 0; i < tableSize; i++)
        {
            item = table + (hash + i) % tableSize;

            // If the element is empty or unused (val == 0) then the item doesn't exist in the table
            if(item->key == -1 || item->val == 0)
            {
                item->key = key;
                item->val = 0;
                break;
            }

            // If it's key matches the input key then return it
            if(item->key == key)
            {
                break;
            }
        }

        // If we landed on a matching key then we're done!
        if(item->key == key)
        {
            return item->val;
        }
        else
        {
            static ElemT outOfMemory;
            return outOfMemory;
        }
    }
};

template<int rows, int cols, int tableSize = cols, class ElemT = float> using SparseMatrix = Matrix<rows, cols, ElemT, Sparse<cols, tableSize, ElemT> >;

//////////////////////////////////////////////////////////// Matrix Minor Memory Delegate ////////////////////////////////////////////////////////////////

template<class ElemT, class MemT> struct Minor
{
    const MemT parent;
    int i, j;

    Minor<ElemT,MemT>(const MemT &obj, int row, int col) : parent(obj), i(row), j(col) { }

    ElemT &operator()(int row, int col) const
    {
        if(row >= i) row++;
        if(col >= j) col++;

        return parent(row,col);
    }
};

//////////////////////////////////////////////////////////// Transpose Delegate ////////////////////////////////////////////////////////////////

template<class ElemT, class MemT> struct Trans
{
    const MemT parent;

    Trans<ElemT,MemT>(const MemT &obj) : parent(obj) { }
    Trans<ElemT,MemT>(const Trans<ElemT,MemT> &obj) : parent(obj.parent) { }

    ElemT &operator()(int row, int col) const
    {
        return parent(col,row);
    }
};

////////////////////////////////////////////////////////// Concatenation Delegates /////////////////////////////////////////////////////////////

template<int leftCols, class ElemT, class LeftMemT, class RightMemT> struct HorzCat
{
    const LeftMemT left;
    const RightMemT right;

    HorzCat<leftCols,ElemT,LeftMemT,RightMemT>(const LeftMemT &l, const RightMemT &r) : left(l), right(r) { }
    HorzCat<leftCols,ElemT,LeftMemT,RightMemT>(const HorzCat<leftCols, ElemT,LeftMemT,RightMemT> &obj) : left(obj.left), right(obj.right) { }

    virtual ~HorzCat<leftCols,ElemT,LeftMemT,RightMemT>() { }

    ElemT &operator()(int row, int col) const
    {
        return col < leftCols? left(row,col) : right(row,col-leftCols);
    }
};

template<int topRows, class ElemT, class TopMemT, class BottomMemT> struct VertCat
{
    const TopMemT top;
    const BottomMemT bottom;

    VertCat<topRows,ElemT,TopMemT,BottomMemT>(const TopMemT &t, const BottomMemT &b) : top(t), bottom(b) { }
    VertCat<topRows,ElemT,TopMemT,BottomMemT>(const VertCat<topRows, ElemT,TopMemT,BottomMemT> &obj) : top(obj.top), bottom(obj.bottom) { }

    virtual ~VertCat<topRows,ElemT,TopMemT,BottomMemT>() { }

    ElemT &operator()(int row, int col) const
    {
        return row < topRows? top(row,col) : bottom(row-topRows,col);
    }
};

#endif // MEMORY_DELEGATE_H
