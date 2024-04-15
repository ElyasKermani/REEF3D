
class lexer;
class chronoWrapper;

class chronoWrapperOuter
{
public:
    chronoWrapperOuter(lexer*);
    ~chronoWrapperOuter();
    void test();
private:
    chronoWrapper* obj;
};