/*
 <DSPHCHRONOLIB>  Copyright (c) 2022 by 
 
 Iván Martínez Estévez (ivan.martinez.estevez@uvigo.es). Universidade de Vigo, Spain.
 Dr José M. Domínguez (jmdominguez@uvigo.es). Universidade de Vigo, Spain.
 Dr Ricardo Canelas (ricardo.canelas@bentley.com). Bentley Systems, Lisbon, Portugal.
 Dr Bonaventura Tagliafierro (btagliafierro@gmail.com). University of Salerno, Italy.
 Dr Orlando García Feal (orlando@uvigo.es). Universidade de Vigo, Spain.
 Professor Alejandro J.C. Crespo (alexbexe@uvigo.es). Universidade de Vigo, Spain.
 Professor Moncho Gómez Gesteira (mggesteira@uvigo.es). Universidade de Vigo, Spain. 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.

 This file is part of DSPHChronoLib. 

 DSPHChronoLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DSPHChronoLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DSPHChronoLib. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JChronoData.cpp \brief Implements the class \ref JChronoData.

#include "JChronoData.h"
#include "Functions.h"
#include "FunDSPHChrono.h"
#include <algorithm>    // std::sort
#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JChValues
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChValues::JChValues(){
  ClassName="JChValues";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JChValues::~JChValues(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JChValues::Reset(){
  LisValue.clear();
}

//==============================================================================
/// Returns data type in text.
//==============================================================================
std::string JChValues::TypeToStr(TpValue type){
  string tx="???";
  switch(type){
    case JChValues::TP_Text:     tx="text";     break;
    case JChValues::TP_Int:      tx="int";      break;
    case JChValues::TP_Uint:     tx="uint";     break;
    case JChValues::TP_Double:   tx="double";   break;
    case JChValues::TP_Int3:     tx="int3";     break;
    case JChValues::TP_Uint3:    tx="uint3";    break;
    case JChValues::TP_Double3:  tx="double3";  break;
  }
  return(tx);
}

//==============================================================================
/// Returns position in LisBody (UINT_MAX when it was not found).
//==============================================================================
unsigned JChValues::IndexByName(const std::string &name)const{
  unsigned ipos=0;
  for(;ipos<GetCount() && LisValue[ipos].name!=name;ipos++);
  return(ipos<GetCount()? ipos: UINT_MAX);
}

//==============================================================================
/// Returns value void.
//==============================================================================
JChValues::StValue JChValues::GetValueVoid(TpValue type,const std::string name)const{
  StValue va;
  va.type=type;
  va.name=name;
  va.vdouble3=TDouble3(0);
  return(va);
}

//==============================================================================
/// Add new value in LisValue.
//==============================================================================
void JChValues::AddValue(const StValue &va){
  if(IndexByName(va.name)!=UINT_MAX)Run_Exceptioon("Name is already in use for other Value.");
  LisValue.push_back(va);
}

//==============================================================================
/// Add value with type string.
//==============================================================================
void JChValues::AddValueStr(const std::string &name,std::string v){
  StValue va=GetValueVoid(TP_Text,name);
  va.vtext=v;
  AddValue(va);
}

//==============================================================================
/// Add value with type int.
//==============================================================================
void JChValues::AddValueInt(const std::string &name,int v){
  StValue va=GetValueVoid(TP_Int,name);
  va.vint=v;
  AddValue(va);
}

//==============================================================================
/// Add value with type unsigned.
//==============================================================================
void JChValues::AddValueUint(const std::string &name,unsigned v){
  StValue va=GetValueVoid(TP_Uint,name);
  va.vuint=v;
  AddValue(va);
}

//==============================================================================
/// Add value with type double.
//==============================================================================
void JChValues::AddValueDouble(const std::string &name,double v){
  StValue va=GetValueVoid(TP_Double,name);
  va.vdouble=v;
  AddValue(va);
}

//==============================================================================
/// Add value with type tint3.
//==============================================================================
void JChValues::AddValueInt3(const std::string &name,tint3 v){
  StValue va=GetValueVoid(TP_Int3,name);
  va.vint3=v;
  AddValue(va);
}

//==============================================================================
/// Add value with type tuint3.
//==============================================================================
void JChValues::AddValueUint3(const std::string &name,tuint3 v){
  StValue va=GetValueVoid(TP_Uint3,name);
  va.vuint3=v;
  AddValue(va);
}

//==============================================================================
/// Add value with type tdouble3.
//==============================================================================
void JChValues::AddValueDouble3(const std::string &name,tdouble3 v){
  StValue va=GetValueVoid(TP_Double3,name);
  va.vdouble3=v;
  AddValue(va);
}

//==============================================================================
/// Checks if the value and type exist.
//==============================================================================
bool JChValues::ExistsValue(const std::string &name,TpValue type)const{
  unsigned ipos=IndexByName(name);
  return(ipos!=UINT_MAX && LisValue[ipos].type==type);
}

//==============================================================================
/// Checks if the value and type exist.
//==============================================================================
bool JChValues::ExistsValue(const std::string &name)const{
  return(IndexByName(name)!=UINT_MAX);
}

//==============================================================================
/// Checks if the value exists and its type is valid. Returns index of value.
//==============================================================================
unsigned JChValues::CheckValueType(const std::string &name,TpValue type,bool optional)const{
  unsigned ipos=IndexByName(name);
  if(ipos==UINT_MAX && !optional)Run_Exceptioon(string("Value with name \'")+name+"\' not found.");
  if(ipos!=UINT_MAX && LisValue[ipos].type!=type)Run_Exceptioon(string("Type of value \'")+name+"\' does not match with type \'"+TypeToStr(type)+"\'.");
  return(ipos);
}

//==============================================================================
/// Returns value with type string.
//==============================================================================
std::string JChValues::GetValueStr(const std::string &name,bool optional,std::string v)const{
  unsigned ipos=CheckValueType(name,TP_Text,optional);
  return(ipos==UINT_MAX?v:LisValue[ipos].vtext);
}

//==============================================================================
/// Returns value with type int.
//==============================================================================
int JChValues::GetValueInt(const std::string &name,bool optional,int v)const{
  unsigned ipos=CheckValueType(name,TP_Int,optional);
  return(ipos==UINT_MAX?v:LisValue[ipos].vint);
}

//==============================================================================
/// Returns value with type unsigned.
//==============================================================================
unsigned JChValues::GetValueUint(const std::string &name,bool optional,unsigned v)const{
  unsigned ipos=CheckValueType(name,TP_Uint,optional);
  return(ipos==UINT_MAX?v:LisValue[ipos].vuint);
}

//==============================================================================
/// Returns value with type double.
//==============================================================================
double JChValues::GetValueDouble(const std::string &name,bool optional,double v)const{
  unsigned ipos=CheckValueType(name,TP_Double,optional);
  return(ipos==UINT_MAX?v:LisValue[ipos].vdouble);
}

//==============================================================================
/// Returns value with type tint3.
//==============================================================================
tint3 JChValues::GetValueInt3(const std::string &name,bool optional,tint3 v)const{
  unsigned ipos=CheckValueType(name,TP_Int3,optional);
  return(ipos==UINT_MAX?v:LisValue[ipos].vint3);
}

//==============================================================================
/// Returns value with type tuint3.
//==============================================================================
tuint3 JChValues::GetValueUint3(const std::string &name,bool optional,tuint3 v)const{
  unsigned ipos=CheckValueType(name,TP_Uint3,optional);
  return(ipos==UINT_MAX?v:LisValue[ipos].vuint3);
}

//==============================================================================
/// Returns value with type tdouble3.
//==============================================================================
tdouble3 JChValues::GetValueDouble3(const std::string &name,bool optional,tdouble3 v)const{
  unsigned ipos=CheckValueType(name,TP_Double3,optional);
  return(ipos==UINT_MAX?v:LisValue[ipos].vdouble3);
}

//==============================================================================
/// Returns value.
//==============================================================================
const JChValues::StValue* JChValues::GetValue(unsigned ipos)const{
  if(ipos>=GetCount())Run_Exceptioon("Number of requested value is invalid.");
  return(&(LisValue[ipos]));
}


//##############################################################################
//# JChBase
//##############################################################################
//==============================================================================
/// Throws exception related to a file or not.
//==============================================================================
void JChBase::RunExceptioon(const std::string &srcfile,int srcline
  ,const std::string &classname,const std::string &method
  ,const std::string &msg,const std::string &file)const
{
  std::string tx;
  if(srcfile.empty())tx=fun::PrintStr("\n*** Exception (%s::%s)\n",ClassName.c_str(),method.c_str());
  else tx=fun::PrintStr("\n*** Exception (%s::%s) at %s:%d\n",ClassName.c_str(),method.c_str(),fun::GetPathLevels(srcfile,3).c_str(),srcline);
  if(!msg.empty())tx=tx+fun::PrintStr("Text: %s\n",msg.c_str());
  if(!file.empty())tx=tx+fun::PrintStr("File: %s\n",file.c_str());
  printf("%s\n",tx.c_str());
  fflush(stdout);
  throw string("#")+tx;
}


//##############################################################################
//# JChBody
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChBody::JChBody(unsigned idb,std::string idname,TpBody type,word mkbound)
  :Idb(idb),IdName(idname),Type(type),MkBound(mkbound)
{
  ClassName="JChBody";
}

//==============================================================================
/// Destructor.
//==============================================================================
JChBody::~JChBody(){
  Reset();
}

//==============================================================================
/// Copy all data from another object.
//==============================================================================
void JChBody::CopyFrom(const JChBody &src){
  ResetRefs();
  Values=src.Values;
  Mass=src.Mass;
  Center=src.Center;
  Inertia=src.Inertia;
  MotionFree=src.MotionFree;
  TranslationFree=src.TranslationFree;
  RotationFree=src.RotationFree;
  LinearVelini=src.LinearVelini;
  AngularVelini=src.AngularVelini;
  ModelFile=src.ModelFile;
  ModelNormal=src.ModelNormal;
  Kfric=src.Kfric;
  Sfric=src.Sfric;
  Restitu=src.Restitu;
  Young=src.Young;
  Poisson=src.Poisson;
  ImposeFric=src.ImposeFric;
  ScaleForce=src.ScaleForce;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JChBody::Reset(){
  ResetRefs();
  Values.Reset();
  Mass=DBL_MAX;
  Center=TDouble3(DBL_MAX);
  Inertia=TMatrix3d(DBL_MAX);
  MotionFree=true;
  TranslationFree=RotationFree=TInt3(1);
  LinearVelini=AngularVelini=TFloat3(0);
  ModelFile="";
  ModelNormal=NorNull;
  Kfric=Sfric=Restitu=Young=Poisson=FLT_MAX;
  ImposeFric=false;
  ScaleForce=TFloat3(FLT_MAX);
}

//==============================================================================
/// Returns data type in text.
//==============================================================================
std::string JChBody::TypeToStr(TpBody type){
  string tx="???";
  switch(type){
    case JChBody::BD_Floating:  tx="Floating";  break;
    case JChBody::BD_Moving:    tx="Moving";    break;
    case JChBody::BD_Fixed:     tx="Fixed";     break;
  }
  return(tx);
}

//==============================================================================
/// Returns model normal in text.
//==============================================================================
std::string JChBody::NormalToStr(TpModelNormal tnor){
  string tx="???";
  switch(tnor){
    case JChBody::NorOriginal:  tx="Original";  break;
    case JChBody::NorInvert:    tx="Invert";    break;
    case JChBody::NorTwoFace:   tx="TwoFace";   break;
  }
  return(tx);
}

//==============================================================================
/// Adds link reference.
//==============================================================================
void JChBody::AddLinkRef(const JChLink* link){
  unsigned ipos=0;
  for(;ipos<GetLinkRefCount() && LinkRefs[ipos]!=link;ipos++);
  if(ipos<GetLinkRefCount())Run_Exceptioon("Number of requested reference is invalid.");
  LinkRefs.push_back((JChLink*)link);
}

//==============================================================================
/// Returns link reference.
//==============================================================================
const JChLink* JChBody::GetLinkRef(unsigned ipos)const{
  if(ipos>=GetLinkRefCount())Run_Exceptioon("Number of requested reference is invalid.");
  return(LinkRefs[ipos]);
}

//==============================================================================
/// Sets parameters for collisions using Simple or Smooth Contacts.
//==============================================================================
void JChBody::SetCollisionData(float kfric,float sfric,float restitu,float young,float poisson){
  Kfric=kfric; Sfric=sfric; Restitu=restitu; Young=young; Poisson=poisson;
}


//==============================================================================
/// Sets initial velocities.
//==============================================================================
void JChBody::SetVelIni(tfloat3 linvelini,tfloat3 angvelini){
  //-Initial velocity.
  LinearVelini=linvelini; 
  AngularVelini=angvelini;
}


//##############################################################################
//# JChBodyFloating
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChBodyFloating::JChBodyFloating(unsigned idb,std::string idname,word mkbound)
  :JChBody(idb,idname,BD_Floating,mkbound)
{
  ClassName="JChBodyFloating";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChBodyFloating::JChBodyFloating(const JChBodyFloating &src)
  :JChBody(src.Idb,src.IdName,BD_Floating,src.MkBound)
{
  ClassName="JChBodyFloating";
  Reset();
  CopyFrom(src);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JChBodyFloating::Reset(){
  JChBody::Reset();
  InputData=false;
  InputFace=InputFomegaAce=TFloat3(0);
  InputLinearVel=InputAngularVel=TFloat3(FLT_MAX);
  OutputCenter=TDouble3(0);
  OutputVel=OutputOmega=TFloat3(0);
}

//==============================================================================
/// Sets floating data.
//==============================================================================
void JChBodyFloating::SetFloatingData(double mass,tdouble3 center,tmatrix3d inertia
  ,tint3 translationfree,tint3 rotationfree,tfloat3 linvelini,tfloat3 angvelini)
{
  Mass=mass;
  Center=center;
  Inertia=inertia;
  //-Motion constraints configuration.
  TranslationFree=translationfree;
  RotationFree=rotationfree;
  MotionFree=(TranslationFree==TInt3(1) && RotationFree==TInt3(1));
  //-Initial velocity. Check if the velocities were initialised before.
  if(LinearVelini ==TFloat3(0))LinearVelini =linvelini; 
  if(AngularVelini==TFloat3(0))AngularVelini=angvelini;
}

//##############################################################################
//# JChBodyMoving
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChBodyMoving::JChBodyMoving(unsigned idb,std::string idname,word mkbound,double mass)
  :JChBody(idb,idname,BD_Moving,mkbound)
{
  ClassName="JChBodyMoving";
  Reset();
  Mass=mass;
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChBodyMoving::JChBodyMoving(const JChBodyMoving &src)
  :JChBody(src.Idb,src.IdName,BD_Moving,src.MkBound)
{
  ClassName="JChBodyMoving";
  Reset();
  CopyFrom(src);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JChBodyMoving::Reset(){
  JChBody::Reset();
  MotionType=MV_None;
  MotionSimple=TDouble3(0);
  MotionMatrix=TMatrix4d(0);
  MotionDt=0;
}

//##############################################################################
//# JChBodyFixed
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChBodyFixed::JChBodyFixed(unsigned idb,std::string idname,word mkbound)
  :JChBody(idb,idname,BD_Fixed,mkbound)
{
  ClassName="JChBodyFixed";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChBodyFixed::JChBodyFixed(const JChBodyFixed &src)
  :JChBody(src.Idb,src.IdName,BD_Fixed,src.MkBound)
{
  ClassName="JChBodyFixed";
  Reset();
  CopyFrom(src);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JChBodyFixed::Reset(){
  JChBody::Reset();
}


//##############################################################################
//# JChLink
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChLink::JChLink(std::string name,TpLink type,unsigned idbody1,unsigned idbody2)
  :Name(name),Type(type),IdBody1(idbody1),IdBody2(idbody2){
  ClassName="JChLink";
}

//==============================================================================
/// Destructor.
//==============================================================================
JChLink::~JChLink(){
  Reset();
}

//==============================================================================
/// Copy all data from another object.
//==============================================================================
void JChLink::CopyFrom(const JChLink &src){
  ResetRefs();
  Values=src.Values;
  Stiffness=src.Stiffness;
  Damping=src.Damping;
  RotPoint=src.RotPoint;
  //RotPoint2=src.RotPoint2;
  RotVector=src.RotVector;
  Pointfb0=src.Pointfb0;
  Pointfb1=src.Pointfb1;
  SlidingVector=src.SlidingVector;
  RotVector2=src.RotVector2;
  RestLength=src.RestLength;
  SvSpring=src.SvSpring;
  CoulombDamping=src.CoulombDamping;
  VariableK=src.VariableK;
  VariableC=src.VariableC;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JChLink::Reset(){
  ResetRefs();
  Values.Reset();
  //RotAngle=DBL_MAX;
  RotPoint=TDouble3(DBL_MAX);
  //RotPoint2=TDouble3(DBL_MAX);
  RotVector=TDouble3(DBL_MAX);
  //Friction=DBL_MAX;
  Stiffness=Damping=0;
  Pointfb0=TDouble3(DBL_MAX);
  Pointfb1=TDouble3(DBL_MAX);
  SlidingVector=TDouble3(DBL_MAX);
  RotVector2=TDouble3(DBL_MAX);
  RestLength=DBL_MAX;
  SvSpring=StrSaveSpring();
  CoulombDamping=0;
  VariableK=VariableC=false;
}

//==============================================================================
/// Returns data type in text.
//==============================================================================
std::string JChLink::TypeToStr(TpLink type){
  string tx="???";
  switch(type){
    case JChLink::LK_Hinge:           tx="Hinge";           break;
    case JChLink::LK_Spheric:         tx="Spheric";         break;
    case JChLink::LK_PointLine:       tx="PointLine";       break;
    case JChLink::LK_LinearSpring:    tx="LinearSpring";    break;
    case JChLink::LK_CoulombDamping:  tx="CoulombDamping";  break;
    case JChLink::LK_Pulley:          tx="Pulley";          break;
    case JChLink::LK_PointFrame:      tx="PointFrame";          break;
  }
  return(tx);
}

//==============================================================================
/// Adds body reference.
//==============================================================================
void JChLink::AddBodyRef(const JChBody* body){
  unsigned ipos=0;
  for(;ipos<GetBodyRefCount() && BodyRefs[ipos]!=body;ipos++);
  if(ipos<GetBodyRefCount())Run_Exceptioon("Number of requested reference is invalid.");
  BodyRefs.push_back((JChBody*)body);
}

//==============================================================================
/// Returns link reference.
//==============================================================================
const JChBody* JChLink::GetBodyRef(unsigned ipos)const{
  if(ipos>=GetBodyRefCount())Run_Exceptioon("Number of requested reference is invalid.");
  return(BodyRefs[ipos]);
}

//##############################################################################
//# JChLinkHinge
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChLinkHinge::JChLinkHinge(std::string name,unsigned idbody1,unsigned idbody2)
  :JChLink(name,LK_Hinge,idbody1,idbody2){
  ClassName="JChLinkHinge";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChLinkHinge::JChLinkHinge(const JChLinkHinge &src)
  :JChLink(src.Name,src.Type,src.IdBody1,src.IdBody2){
  ClassName="JChLinkHinge";
  Reset();
  CopyFrom(src);
}

//##############################################################################
//# JChLinkSpheric
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChLinkSpheric::JChLinkSpheric(std::string name,unsigned idbody1,unsigned idbody2)
  :JChLink(name,LK_Spheric,idbody1,idbody2){
  ClassName="JChLinkSpheric";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChLinkSpheric::JChLinkSpheric(const JChLinkSpheric &src)
  :JChLink(src.Name,src.Type,src.IdBody1,src.IdBody2){
  ClassName="JChLinkSpheric";
  Reset();
  CopyFrom(src);
}

//##############################################################################
//# JChLinkPointLine
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChLinkPointLine::JChLinkPointLine(std::string name,unsigned idbody1,unsigned idbody2)
  :JChLink(name,LK_PointLine,idbody1,idbody2){
  ClassName="JChLinkPointLine";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChLinkPointLine::JChLinkPointLine(const JChLinkPointLine &src)
  :JChLink(src.Name,src.Type,src.IdBody1,src.IdBody2){
  ClassName="JChLinkPointLine";
  Reset();
  CopyFrom(src);
}

//##############################################################################
//# JChLinkLinearSpring
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChLinkLinearSpring::JChLinkLinearSpring(std::string name,unsigned idbody1,unsigned idbody2)
  :JChLink(name,LK_LinearSpring,idbody1,idbody2){
  ClassName="JChLinkLinearSpring";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChLinkLinearSpring::JChLinkLinearSpring(const JChLinkLinearSpring &src)
  :JChLink(src.Name,src.Type,src.IdBody1,src.IdBody2){
  ClassName="JChLinkLinearSpring";
  Reset();
  CopyFrom(src);
}

//##############################################################################
//# JChLinkCoulombDamping
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChLinkCoulombDamping::JChLinkCoulombDamping(std::string name,unsigned idbody1,unsigned idbody2)
  :JChLink(name,LK_CoulombDamping,idbody1,idbody2){
  ClassName="JChLinkCoulombDamping";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChLinkCoulombDamping::JChLinkCoulombDamping(const JChLinkCoulombDamping &src)
  :JChLink(src.Name,src.Type,src.IdBody1,src.IdBody2){
  ClassName="JChLinkCoulombDamping";
  Reset();
  CopyFrom(src);
}

//##############################################################################
//# JChLinkPulley
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChLinkPulley::JChLinkPulley(std::string name,unsigned idbody1,unsigned idbody2)
  :JChLink(name,LK_Pulley,idbody1,idbody2){
  ClassName="JChLinkPulley";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChLinkPulley::JChLinkPulley(const JChLinkPulley &src)
  :JChLink(src.Name,src.Type,src.IdBody1,src.IdBody2){
  ClassName="JChLinkPulley";
  Reset();
  CopyFrom(src);
}
//==============================================================================
/// Copy all data from another object.
//==============================================================================
void JChLinkPulley::CopyFrom(const JChLinkPulley &src){
  ResetRefs();
  Radius=src.Radius;
  Radius2=src.Radius2;
  JChLink::CopyFrom(src);
}

//##############################################################################
//# JChronoData
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChronoData::JChronoData(){
  ClassName="JChronoData";
  Reset();
}

//==============================================================================
/// Constructor for copies.
//==============================================================================
JChronoData::JChronoData(const JChronoData &src){
  ClassName="JChronoData";
  Reset();
  *this=src;
}

//==============================================================================
/// Destructor.
//==============================================================================
JChronoData::~JChronoData(){
  Reset();
}

//==============================================================================
/// Overload assignment operator.
//==============================================================================
JChronoData& JChronoData::operator=(const JChronoData &src){
  if(this!=&src){
    Reset();
    //-Copies general data.
    DataDir=src.DataDir;
    UseNSCChrono=src.UseNSCChrono;
    Dp=src.Dp;
    CollisionDp=src.CollisionDp;
    OmpThreads=src.OmpThreads;
    Solver=src.Solver;
    TimeStepper=src.TimeStepper;
    ExecMode=src.ExecMode;
    ContactMethod=src.ContactMethod;
    UseVariableCoeff=src.UseVariableCoeff;
    UseCollision=src.UseCollision;
    UseGravity=src.UseGravity;
    Gravity=src.Gravity;
    FtPause=src.FtPause;
    //-Copies bodies.
    for(unsigned c=0;c<src.GetBodyCount();c++){
      const JChBody *body=src.LisBody[c];
      if(body->Type==JChBody::BD_Floating){
        JChBodyFloating *newbody=new JChBodyFloating(*(JChBodyFloating*)body);
        LisBody.push_back(newbody);
      }
      else if(body->Type==JChBody::BD_Moving){
        JChBodyMoving *newbody=new JChBodyMoving(*(JChBodyMoving*)body);
        LisBody.push_back(newbody);
      }
      else if(body->Type==JChBody::BD_Fixed){
        JChBodyFixed *newbody=new JChBodyFixed(*(JChBodyFixed*)body);
        LisBody.push_back(newbody);
      }
      else Run_Exceptioon("Class invalid for data copying body.");
    }
    //-Copies links.
    for(unsigned c=0;c<src.GetLinkCount();c++){
      const JChLink *link=src.LisLink[c];
      if(link->Type==JChLink::LK_Hinge){
        JChLinkHinge *newlink=new JChLinkHinge(*(JChLinkHinge*)link);
        LisLink.push_back(newlink);
      }
      else if(link->Type==JChLink::LK_Pulley){
        JChLinkPulley *newlink=new JChLinkPulley(*(JChLinkPulley*)link);
        LisLink.push_back(newlink);
      }
      else if(link->Type==JChLink::LK_Spheric){
        JChLinkSpheric *newlink=new JChLinkSpheric(*(JChLinkSpheric*)link);
        LisLink.push_back(newlink);
      }
      else if(link->Type==JChLink::LK_PointLine){
        JChLinkPointLine *newlink=new JChLinkPointLine(*(JChLinkPointLine*)link);
        LisLink.push_back(newlink);
      }
      else if(link->Type==JChLink::LK_LinearSpring){
        JChLinkLinearSpring *newlink=new JChLinkLinearSpring(*(JChLinkLinearSpring*)link);
        LisLink.push_back(newlink);
      }
      else if(link->Type==JChLink::LK_CoulombDamping){
        JChLinkCoulombDamping *newlink=new JChLinkCoulombDamping(*(JChLinkCoulombDamping*)link);
        LisLink.push_back(newlink);
      }
      else Run_Exceptioon("Class invalid for copying link.");
    }

    Prepare();
  }
  return(*this);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JChronoData::Reset(){
  DataDir="";
  UseNSCChrono=false;
  Dp=0;
  CollisionDp=0.5;
  for(unsigned c=0;c<LisBody.size();c++)delete LisBody[c];
  for(unsigned c=0;c<LisLink.size();c++)delete LisLink[c];
  LisBody.clear();
  LisLink.clear();
  Solver=0;
  TimeStepper=0;
  OmpThreads=1;
  ExecMode="";
  ContactMethod=NSC;
  UseVariableCoeff=false;
  UseCollision=false;
  UseGravity=false;
  Gravity=TFloat3(0);
  FtPause=0;
}

//==============================================================================
/// Prepares and checks references.
//==============================================================================
void JChronoData::Prepare(){
  //-Configures bodies with its references to links.
  for(unsigned cb=0;cb<LisBody.size();cb++){
    JChBody *body=LisBody[cb];
    body->ResetRefs();
    for(unsigned ck=0;ck<LisLink.size();ck++){
      if(LisLink[ck]->IdBody1==body->Idb || LisLink[ck]->IdBody2==body->Idb)body->AddLinkRef(LisLink[ck]);
    }
  }
  //-Configures links with its references to bodies.
  for(unsigned ck=0;ck<LisLink.size();ck++){
    JChLink *link=LisLink[ck];
    link->ResetRefs();
    unsigned ipos1=BodyIndexById(link->IdBody1);
    unsigned ipos2=BodyIndexById(link->IdBody2);
    if(ipos1!=UINT_MAX)link->AddBodyRef(LisBody[ipos1]);
    else Run_Exceptioon("Link has reference to body invalid.");
    if(ipos2!=UINT_MAX)link->AddBodyRef(LisBody[ipos2]);
  }
}

//==============================================================================
/// Prepares and checks references.
//==============================================================================
void JChronoData::CheckData(){
  if(!LisBody.size())Run_Exceptioon("There are no defined bodies to use Chrono.");
  //-Checks floating body data.
  for(unsigned cb=0;cb<LisBody.size();cb++){
    if(LisBody[cb]->Type==JChBody::BD_Floating){
      const JChBodyFloating *body=(const JChBodyFloating *)LisBody[cb];
      if(body->GetMass()==DBL_MAX || body->GetCenter()==TDouble3(DBL_MAX) || body->GetInertia()==TMatrix3d(DBL_MAX))
        Run_Exceptioon(string("BodyFloating \'")+body->IdName+"\' without data from floating body with mkbound.");
    }
    else if(LisBody[cb]->Type==JChBody::BD_Moving){
      const JChBodyMoving *body=(const JChBodyMoving *)LisBody[cb];
      if(body->GetMass()==DBL_MAX || body->GetCenter()==TDouble3(DBL_MAX))
        Run_Exceptioon(string("BodyMoving \'")+body->IdName+"\' without mass or center data.");
    }
  }
}

//==============================================================================
/// Returns error for adding a new body.
//==============================================================================
std::string JChronoData::CheckAddBodyError(unsigned idb,std::string idname,word mkbound)const{
  if(idname=="NULL")return(fun::PrintStr("Id=\'%s\' is invalid.",idname.c_str()));
  if(BodyIndexById(idb)!=UINT_MAX)return(fun::PrintStr("Idb code %u is already in use for other body.",idb));
  if(BodyIndexByName(idname)!=UINT_MAX)return(fun::PrintStr("Id=\'%s\' is already in use for other body.",idname.c_str()));
  if(BodyIndexByMkBound(mkbound)!=UINT_MAX)return(fun::PrintStr("MkBound %u is already in use for other body with id=\'%s\'.",mkbound,LisBody[BodyIndexByMkBound(mkbound)]->IdName.c_str()));
  return("");
}

//==============================================================================
/// Add new BodyFloating.
//==============================================================================
JChBodyFloating* JChronoData::AddBodyFloating(unsigned idb,std::string idname,word mkbound,std::string fileinfo){
  JChBodyFloating* body=NULL;
  const string errtx=CheckAddBodyError(idb,idname,mkbound);
  if(!errtx.empty())Run_ExceptioonFile(errtx,fileinfo);
  body=new JChBodyFloating(idb,idname,mkbound);
  LisBody.push_back(body);
  return(body);
}

//==============================================================================
/// Add new BodyMoving.
//==============================================================================
JChBodyMoving* JChronoData::AddBodyMoving(unsigned idb,std::string idname,word mkbound,double mass,std::string fileinfo){
  JChBodyMoving* body=NULL;
  const string errtx=CheckAddBodyError(idb,idname,mkbound);
  if(!errtx.empty())Run_ExceptioonFile(errtx,fileinfo);
  body=new JChBodyMoving(idb,idname,mkbound,mass);
  LisBody.push_back(body);
  return(body);
}

//==============================================================================
/// Add new BodyFixed.
//==============================================================================
JChBodyFixed* JChronoData::AddBodyFixed(unsigned idb,std::string idname,word mkbound,std::string fileinfo){
  JChBodyFixed* body=NULL;
  const string errtx=CheckAddBodyError(idb,idname,mkbound);
  if(!errtx.empty())Run_ExceptioonFile(errtx,fileinfo);
  body=new JChBodyFixed(idb,idname,mkbound);
  LisBody.push_back(body);
  return(body);
}

//==============================================================================
/// Add new LinkHinge.
//==============================================================================
JChLinkHinge* JChronoData::AddLinkHinge(std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo){
  JChLinkHinge* link=NULL;
  if(LinkIndexByName(name)!=UINT_MAX)Run_ExceptioonFile(fun::PrintStr("Name=\'%s\' is already in use for other Link.",name.c_str()),fileinfo);
  link=new JChLinkHinge(name,idbody1,idbody2);
  LisLink.push_back(link);
  return(link);
}

//==============================================================================
/// Add new LinkPulley.
//==============================================================================
JChLinkPulley* JChronoData::AddLinkPulley(std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo){
  JChLinkPulley* link=NULL;
  if(LinkIndexByName(name)!=UINT_MAX)Run_ExceptioonFile(fun::PrintStr("Name=\'%s\' is already in use for other Link.",name.c_str()),fileinfo);
  link=new JChLinkPulley(name,idbody1,idbody2);
  LisLink.push_back(link);
  return(link);
}

//==============================================================================
/// Add new LinkSpheric.
//==============================================================================
JChLinkSpheric* JChronoData::AddLinkSpheric(std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo){
  JChLinkSpheric* link=NULL;
  if(LinkIndexByName(name)!=UINT_MAX)Run_ExceptioonFile(fun::PrintStr("Name=\'%s\' is already in use for other Link.",name.c_str()),fileinfo);
  link=new JChLinkSpheric(name,idbody1,idbody2);
  LisLink.push_back(link);
  return(link);
}

//==============================================================================
/// Add new LinkPointLine.
//==============================================================================
JChLinkPointLine* JChronoData::AddLinkPointLine(std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo){
  JChLinkPointLine* link=NULL;
  if(LinkIndexByName(name)!=UINT_MAX)Run_ExceptioonFile(fun::PrintStr("Name=\'%s\' is already in use for other Link.",name.c_str()),fileinfo);
  link=new JChLinkPointLine(name,idbody1,idbody2);
  LisLink.push_back(link);
  return(link);
}

//==============================================================================
/// Add new LinkLinearSpring.
//==============================================================================
JChLinkLinearSpring* JChronoData::AddLinkLinearSpring(std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo){
  JChLinkLinearSpring* link=NULL;
  if(LinkIndexByName(name)!=UINT_MAX)Run_ExceptioonFile(fun::PrintStr("Name=\'%s\' is already in use for other Link.",name.c_str()),fileinfo);
  link=new JChLinkLinearSpring(name,idbody1,idbody2);
  LisLink.push_back(link);
  return(link);
}

//==============================================================================
/// Add new LinkCoulombDamping.
//==============================================================================
JChLinkCoulombDamping* JChronoData::AddLinkCoulombDamping(std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo){
  JChLinkCoulombDamping* link=NULL;
  if(LinkIndexByName(name)!=UINT_MAX)Run_ExceptioonFile(fun::PrintStr("Name=\'%s\' is already in use for other Link.",name.c_str()),fileinfo);
  link=new JChLinkCoulombDamping(name,idbody1,idbody2);
  LisLink.push_back(link);
  return(link);
}

//==============================================================================
/// Returns position in LisBody (UINT_MAX when it was not found).
//==============================================================================
unsigned JChronoData::BodyIndexById(unsigned idb)const{
  unsigned ipos=0;
  for(;ipos<GetBodyCount()&&LisBody[ipos]->Idb!=idb;ipos++);
  return(ipos<GetBodyCount()?ipos:UINT_MAX);
}

//==============================================================================
/// Returns position in LisBody (UINT_MAX when it was not found).
//==============================================================================
unsigned JChronoData::BodyIndexByMkBound(word mkbound)const{
  unsigned ipos=0;
  for(;ipos<GetBodyCount()&&LisBody[ipos]->MkBound!=mkbound;ipos++);
  return(ipos<GetBodyCount()?ipos:UINT_MAX);
}

//==============================================================================
/// Returns position in LisBody (UINT_MAX when it was not found).
//==============================================================================
unsigned JChronoData::BodyIndexByName(const std::string &name)const{
  unsigned ipos=0;
  for(;ipos<GetBodyCount()&&LisBody[ipos]->IdName!=name;ipos++);
  return(ipos<GetBodyCount()?ipos:UINT_MAX);
}

//==============================================================================
/// Returns position in LisLink (UINT_MAX when it was not found).
//==============================================================================
unsigned JChronoData::LinkIndexByName(const std::string &name)const{
  unsigned ipos=0;
  for(;ipos<GetLinkCount()&&LisLink[ipos]->Name!=name;ipos++);
  return(ipos<GetLinkCount()?ipos:UINT_MAX);
}

//==============================================================================
/// Returns the requested body.
//==============================================================================
const JChBody* JChronoData::GetBody(unsigned ipos)const{
  if(ipos>=GetBodyCount())Run_Exceptioon("Number of requested body is invalid.");
  return(LisBody[ipos]);
}

//==============================================================================
/// Returns the requested body by mkbound.
//==============================================================================
const JChBody* JChronoData::GetBodyByMk(word mkbound)const{
  const unsigned ipos=BodyIndexByMkBound(mkbound);
  if(ipos>=GetBodyCount())Run_Exceptioon("Mk of requested body is invalid.");
  return(LisBody[ipos]);
}

//==============================================================================
/// Returns true if the requested body belongs to any TpLink
//==============================================================================
bool JChronoData::BodyBelongsLink(const unsigned idBody1,JChLink::TpLink type)const{
  std::vector<const JChLink*> links_tp=GetLinkByTp(type);
  for(unsigned l=0;l<links_tp.size();l++){
    const JChLink* link=links_tp[l];
    if(idBody1==link->IdBody1 || idBody1==link->IdBody2)return(true);
  }
  return(false);
}

//==============================================================================
/// Returns the requested link.
//==============================================================================
std::vector<const JChLink*> JChronoData::GetLinkByTp(JChLink::TpLink type)const{
  std::vector<const JChLink*> links_tp;
  for(unsigned l=0;l<GetLinkCount();l++)if(GetLink(l)->Type==type)links_tp.push_back(GetLink(l));
  return(links_tp);
}

//==============================================================================
/// Returns the requested link.
//==============================================================================
JChLink* JChronoData::GetLink(unsigned ipos)const{
  if(ipos>=GetLinkCount())Run_Exceptioon("Number of requested link is invalid.");
  return(LisLink[ipos]);
}

//==============================================================================
/// Returns the position of the requested link.
//==============================================================================
unsigned JChronoData::GetPosLink(const JChLink *link)const{
  if(!link)Run_Exceptioon("The link object is missing.");
  unsigned l;
  for(l=0;l<GetLinkCount() && LisLink[l]!=link;l++);
  if(l>=GetLinkCount())Run_Exceptioon("The requested link was not found.");
  return(l);
}

//==============================================================================
/// Returns bodyfloating object with indicated mkbound.
//==============================================================================
const JChBodyFloating* JChronoData::GetBodyFloating(word mkbound)const{
  const unsigned ipos=BodyIndexByMkBound(mkbound);
  if(ipos<GetBodyCount() && LisBody[ipos]->Type!=JChBody::BD_Floating)Run_Exceptioon("Body with indicated Mk is not a floating body.");
  return(ipos==UINT_MAX?NULL:(const JChBodyFloating*)LisBody[ipos]);
}

//==============================================================================
/// Returns bodymoving object with indicated mkbound.
//==============================================================================
const JChBodyMoving* JChronoData::GetBodyMoving(word mkbound)const{
  const unsigned ipos=BodyIndexByMkBound(mkbound);
  if(ipos<GetBodyCount() && LisBody[ipos]->Type!=JChBody::BD_Moving)Run_Exceptioon("Body with indicated Mk is not a moving body.");
  return(ipos==UINT_MAX?NULL:(const JChBodyMoving*)LisBody[ipos]);
}

//==============================================================================
/// Returns bodyfixed object with indicated mkbound.
//==============================================================================
const JChBodyFixed* JChronoData::GetBodyFixed(word mkbound)const{
  const unsigned ipos=BodyIndexByMkBound(mkbound);
  if(ipos<GetBodyCount() && LisBody[ipos]->Type!=JChBody::BD_Fixed)Run_Exceptioon("Body with indicated Mk is not a fixed body.");
  return(ipos==UINT_MAX?NULL:(const JChBodyFixed*)LisBody[ipos]);
}

////==============================================================================
/////Returns the required node
////==============================================================================
//const JChronoData::StNode* JChronoData::GetStNode(const unsigned n)const{
//  if(!LisNode.size())fun::Run_ExceptioonFun("Cannot access to the requiered using the \'Automate\' mode. You should create manually a set of <nodes> instead of using <pointA> and <pointB>");
//  if(n<0||n>=(unsigned)LisNode.size())fun::Run_ExceptioonFun(fun::PrintStr("Index node is out of range. Expected: [0,%d]. Sent: %d.",(unsigned)LisNode.size()-1,n));
//  return(LisNode[n]);
//}

//==============================================================================
///Converts the enum type to string to show the solver
//==============================================================================
std::string JChronoData::SolverToStr()const{
  std::string tx="???";
  switch (Solver){
    case TpSolverType::BB:       tx="Barzilai-Borwein";  break;
    case TpSolverType::MINRES:   tx="MINimum RESidual";  break;
  }
  return (tx);
}

//==============================================================================
///Converts the enum type to string to show the timestepper
//==============================================================================
std::string JChronoData::TimeStepperToStr()const{
  std::string tx="???";
  switch (TimeStepper){
    case TpTStepperType::EULER_IL:  tx="Euler Implicit Linearized";   break;
  }
  return (tx);
}

//==============================================================================
///Converts the contact method type to string
//==============================================================================
std::string JChronoData::ContactMethodToStr()const{
  std::string tx="???";
  switch (ContactMethod){
    case TpContactMethod::NSC:  tx="Non-Smooth Contacts (NSC)";  break;
    case TpContactMethod::SMC:  tx="SMooth Contacts (SMC)";  break;
  }
  return (tx);
}