#include <memory>
#include <set>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/SharedTracks.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "RecoVertex/MultiVertexFit/interface/MultiVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

template <class VTX>
class TemplatedVertexMerger : public edm::stream::EDProducer<> {
    public:
	typedef std::vector<VTX> Product;
	TemplatedVertexMerger(const edm::ParameterSet &params);

	virtual void produce(edm::Event &event, const edm::EventSetup &es) override;

    private:
       bool trackFilter(const reco::TrackRef &track) const;
  
       void svMerger(typename Product::iterator sv, typename Product::iterator sv2);
       void checkMergedVtx(typename Product::iterator svtx, edm::ESHandle<TransientTrackBuilder> theTTBuilder, const reco::Vertex pv);

       edm::EDGetTokenT<reco::VertexCollection> token_primaryVertex;
       edm::EDGetTokenT<Product> 	        token_secondaryVertex;
       double					maxFraction;
       double					minSignificance;
       bool                                     doMerging;
};

template <class VTX>
TemplatedVertexMerger<VTX>::TemplatedVertexMerger(const edm::ParameterSet &params) :
	maxFraction(params.getParameter<double>("maxFraction")),
	minSignificance(params.getParameter<double>("minSignificance")),
	doMerging(params.getParameter<bool>("doMerging"))
{
        token_primaryVertex = consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("primaryVertices"));
	token_secondaryVertex = consumes<Product>(params.getParameter<edm::InputTag>("secondaryVertices"));
	produces<Product>();
}

template <class VTX>
void TemplatedVertexMerger<VTX>::produce(edm::Event &event, const edm::EventSetup &es)
{
	using namespace reco;

	edm::Handle<VertexCollection> primaryVertices;
	event.getByToken(token_primaryVertex, primaryVertices);

	edm::Handle<Product> secondaryVertices;
	event.getByToken(token_secondaryVertex, secondaryVertices);

	edm::ESHandle<TransientTrackBuilder> theTTBuilder;
	es.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
 
	const reco::Vertex &pv = (*primaryVertices)[0];
	
        VertexDistance3D dist;
	std::auto_ptr<Product> recoVertices(new Product);
	for(typename Product::const_iterator sv = secondaryVertices->begin(); 
	    sv != secondaryVertices->end(); ++sv) {
	  recoVertices->push_back(*sv);
        }
	
	
	if(doMerging){ // Allowing merging instead of dropping one of the secondary vertices.
          
	  for(typename Product::iterator sv = recoVertices->begin();
	      sv < recoVertices->end(); ++sv)
	    {
	      VertexState s1(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error()));
	      for(typename Product::iterator sv2 = sv; sv2 < recoVertices->end(); ++sv2) {
		
		VertexState s2(RecoVertex::convertPos(sv2->position()),RecoVertex::convertError(sv2->error()));
		double fr=vertexTools::computeSharedTracks(*sv, *sv2);
		double fr2=vertexTools::computeSharedTracks(*sv2, *sv);
		
		if (fr > maxFraction && dist.distance(s1,s2).significance() < minSignificance && fr >= fr2 && sv-sv2!=0)
		  {
		    svMerger(sv,sv2); // Merging sv2 with sv and removing sv2
		    checkMergedVtx(sv,theTTBuilder,pv); // Testing if the new vertex make sense. Print out 
		    sv2=recoVertices->erase(sv2)-1;
		  }
		else if (fr2 > maxFraction && dist.distance(s1,s2).significance() < minSignificance && sv-sv2!=0)
		  {
		    svMerger(sv2,sv);  // Merging sv wih sv2 and removing sv
		    checkMergedVtx(sv2,theTTBuilder,pv); // Testing if the new vertex make sense. Print out              
		    sv=recoVertices->erase(sv)-1;
		    sv2=recoVertices->end();
		  }
	      }
	    }
	}
	
	else { 	// Just dropping sv if sv and sv2 shared
          for(typename Product::iterator sv = recoVertices->begin();
	      sv != recoVertices->end(); ++sv) 
	    {
	      
	      bool shared=false;
	      VertexState s1(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error()));
	      for(typename Product::iterator sv2 = recoVertices->begin();
		  sv2 != recoVertices->end(); ++sv2) {
		VertexState s2(RecoVertex::convertPos(sv2->position()),RecoVertex::convertError(sv2->error()));
		double fr=vertexTools::computeSharedTracks(*sv2, *sv);
		if (fr > maxFraction && dist.distance(s1,s2).significance() < minSignificance && sv-sv2!=0 
                    && fr >= vertexTools::computeSharedTracks(*sv, *sv2) )
 		  {
		    shared=true; 
		  }
	      }
              if(shared) { sv=recoVertices->erase(sv)-1; }
	    }
	}
	
	event.put(recoVertices);
	
} //End produce

// Check the merged vertex
template <>
void TemplatedVertexMerger<reco::Vertex>::checkMergedVtx(typename Product::iterator svtx, edm::ESHandle<TransientTrackBuilder> theTTBuilder, const reco::Vertex pv)
{
  VertexDistance3D vdist;
  VertexDistanceXY vdist2d;
  
  AdaptiveVertexFitter theAdaptiveFitter(
					 //GeometricAnnealing(fitterSigmacut, fitterTini, fitterRatio),
					 GeometricAnnealing(3, 256, 0.25),
					 DefaultLinearizationPointFinder(),
					 KalmanVertexUpdator<5>(),
					 KalmanVertexTrackCompatibilityEstimator<5>(),
					 KalmanVertexSmoother() );
  
  std::vector<reco::TransientTrack> tt_vtx;
  
  for(reco::Vertex::trackRef_iterator t_svtx = svtx->tracks_begin(); t_svtx!= svtx->tracks_end(); ++t_svtx){
    reco::Track rftrk = svtx->refittedTrack(*t_svtx);
    reco::TransientTrack transientTrack = theTTBuilder->build(rftrk);
    tt_vtx.push_back(transientTrack);
  }

  TransientVertex testFitVertex;
  testFitVertex = theAdaptiveFitter.vertex(tt_vtx);
  
  if(testFitVertex.isValid()){
    
    //std::cout << "Original tracks merged vtx: " << svtx->tracksSize() << "  Refitted tracks merged vtx: " <<  testFitVertex.refittedTracks().size() <<  std::endl;
    float diffsize = fabs(svtx->tracksSize() - testFitVertex.refittedTracks().size());
    
    Measurement1D dlen= vdist.distance(pv,testFitVertex);
    Measurement1D dlen2= vdist2d.distance(pv,testFitVertex);
    
    GlobalVector dir;
    std::vector<reco::TransientTrack> ts = testFitVertex.originalTracks();
    for(std::vector<reco::TransientTrack>::const_iterator i = ts.begin();
	i != ts.end(); ++i) {
      float w = testFitVertex.trackWeight(*i);
      if (w > 0.5) dir+=i->impactPointState().globalDirection();
    }
 
  
    GlobalPoint psv(testFitVertex.position().x(),testFitVertex.position().y(),testFitVertex.position().z());
    GlobalPoint ppv(pv.position().x(),pv.position().y(),pv.position().z());
    
    float vscal = dir.unit().dot((psv-ppv).unit());
 
    //std::cout <<"-------------------------------------------------------------------------------------------------------------------------------------"<<  std::endl;
    //std::cout << "dlen: " << dlen.value() << " dlen2: " << dlen2.value()   <<  std::endl;
    //std::cout << "dlen Sig: " << dlen.significance() << " dlen2 Sig: " << dlen2.significance() <<" normChiSq sv:  " << testFitVertex.normalisedChiSquared() << " vscal: " << vscal  <<  std::endl;
    //std::cout <<testFitVertex.position().x()<<" "<<testFitVertex.position().y()<<" "<<testFitVertex.position().z()<<std::endl;
    //std::cout <<pv.position().x()<<" "<<pv.position().y()<<" "<<pv.position().z()<<std::endl;
    
  
    if(!(diffsize < 2 && dlen.significance() > 0.5 && vscal > 0.95 && testFitVertex.normalisedChiSquared() < 10 && dlen2.significance() > 2.5)){
      //std::cout << "WARNING: Merged vertex not converging or not displaced "<<  std::endl;
      //std::cout <<"-------------------------------------------------------------------------------------------------------------------------------------"<<  std::endl; 
      //std::cout << "dlen: " << dlen.value() << " dlen2: " << dlen2.value()   <<  std::endl;                                                              
      //std::cout << "dlen Sig: " << dlen.significance() << " dlen2 Sig: " << dlen2.significance() <<" normChiSq sv:  " << testFitVertex.normalisedChiSquared() << " vscal: " << vscal <<  std::endl;                                                                                                                                                               
      //std::cout <<"Diff Size:  "<<diffsize<<std::endl;                                                       
    }
  }else std::cout << "WARNING: Merged vertex not valid"<<  std::endl;  

}
//------
template <>
void TemplatedVertexMerger<reco::VertexCompositePtrCandidate>::checkMergedVtx(typename Product::iterator svtx, edm::ESHandle<TransientTrackBuilder> theTTBuilder, const reco::Vertex pv)
{
  
  VertexDistance3D vdist;
  VertexDistanceXY vdist2d;
  
  AdaptiveVertexFitter theAdaptiveFitter(
					 //GeometricAnnealing(fitterSigmacut, fitterTini, fitterRatio),
					 GeometricAnnealing(3, 256, 0.25),
					 DefaultLinearizationPointFinder(),
					 KalmanVertexUpdator<5>(),
					 KalmanVertexTrackCompatibilityEstimator<5>(),
					 KalmanVertexSmoother() );
  
  std::vector<reco::TransientTrack> tt_vtx;
   for(std::vector<reco::CandidatePtr>::const_iterator t_svtx = svtx->daughterPtrVector().begin(); t_svtx != svtx->daughterPtrVector().end(); ++t_svtx){
     reco::TransientTrack transientTrack = theTTBuilder->build(*t_svtx);
     tt_vtx.push_back(transientTrack);
  }
  
  TransientVertex testFitVertex;
  testFitVertex = theAdaptiveFitter.vertex(tt_vtx);
  
  if(testFitVertex.isValid()){
    float diffsize = fabs(svtx->daughterPtrVector().size()- testFitVertex.refittedTracks().size());
    //std::cout << "Original tracks merged vtx: " << svtx->daughterPtrVector().size() << "  Refitted tracks merged vtx: " <<  testFitVertex.refittedTracks().size()  <<  std::endl;
  
    Measurement1D dlen= vdist.distance(pv,testFitVertex);
    Measurement1D dlen2= vdist2d.distance(pv,testFitVertex);
  
    GlobalVector dir;
    std::vector<reco::TransientTrack> ts = testFitVertex.originalTracks();
    for(std::vector<reco::TransientTrack>::const_iterator i = ts.begin();
	i != ts.end(); ++i) {
      float w = testFitVertex.trackWeight(*i);
      if (w > 0.5) dir+=i->impactPointState().globalDirection();
    }
    
    GlobalPoint psv(testFitVertex.position().x(),testFitVertex.position().y(),testFitVertex.position().z());
    GlobalPoint ppv(pv.position().x(),pv.position().y(),pv.position().z());
    float vscal = dir.unit().dot((psv-ppv).unit());

     if(!(diffsize < 2 && dlen.significance() > 0.5 && vscal > 0.95 && testFitVertex.normalisedChiSquared() < 10 && dlen2.significance() > 2.5)){
      std::cout << "WARNING: Merged vertex not converging or not displaced "<<  std::endl;
    }
  
    //std::cout << "dlen: " << dlen.value() << " dlen2: " << dlen2.value()   <<  std::endl;
    //std::cout << "dlen Sig: " << dlen.significance() << " dlen2 Sig: " << dlen2.significance() <<" normChiSq sv:  " << testFitVertex.normalisedChiSquared() << " vscal: " << vscal  <<std::endl;

  }else std::cout << "WARNING: Merged vertex not valid"<<  std::endl;
  
}

// Merger. Vertex svi becomes the merging of svi and svj ---------- 
template <>
void TemplatedVertexMerger<reco::Vertex>::svMerger(typename Product::iterator svi, typename Product::iterator svj)
{
  bool sharingTracks = false;
  
  for(reco::Vertex::trackRef_iterator ti_sv = svj->tracks_begin(); ti_sv!= svj->tracks_end(); ++ti_sv){
    reco::Vertex::trackRef_iterator it = find(svi->tracks_begin(), svi->tracks_end(), *ti_sv);
    if (it==svi->tracks_end()){
      svi->add( *ti_sv, svj->refittedTrack(*ti_sv), svj->trackWeight(*ti_sv));
    }else sharingTracks=true;
  }
  
  if(sharingTracks){
    // create backup track containers from the main vertex
    std::vector<reco::TrackBaseRef> tracks_;
    std::vector<reco::Track> refittedTracks_;
    std::vector<float> weights_;
    
    for(reco::Vertex::trackRef_iterator it = svi->tracks_begin(); it!=svi->tracks_end(); it++) {
      tracks_.push_back( *it);
      refittedTracks_.push_back( svi->refittedTrack(*it));
      weights_.push_back( svi->trackWeight(*it) );
    }
    
    // delete tracks and add all tracks back, and check in which vertex the weight is larger
    svi->removeTracks();
    std::vector<reco::Track>::iterator it2 = refittedTracks_.begin();
    std::vector<float>::iterator it3 = weights_.begin();
    for(reco::Vertex::trackRef_iterator it = tracks_.begin(); it!=tracks_.end(); it++, it2++, it3++){
      float weight = *it3;
      float weight2= svj->trackWeight(*it);
      reco::Track refittedTrackWithLargerWeight = *it2;
      if( weight2 >weight) {
	weight = weight2;
	refittedTrackWithLargerWeight = svj->refittedTrack(*it);
      }
      svi->add(*it , refittedTrackWithLargerWeight , weight);
    }
  }
  
}

template <>
void TemplatedVertexMerger<reco::VertexCompositePtrCandidate>::svMerger(typename Product::iterator svi, typename Product::iterator svj)
{
  for(std::vector<reco::CandidatePtr>::const_iterator ti_sv = svj->daughterPtrVector().begin(); ti_sv != svj->daughterPtrVector().end(); ++ti_sv){
    std::vector<reco::CandidatePtr>::const_iterator it = find(svi->daughterPtrVector().begin(), svi->daughterPtrVector().end(), *ti_sv);
    if (it==svi->daughterPtrVector().end()){ 
      svi->addDaughter(*ti_sv);
    }
  }
}


typedef TemplatedVertexMerger<reco::Vertex> VertexMerger;
typedef TemplatedVertexMerger<reco::VertexCompositePtrCandidate> CandidateVertexMerger;

DEFINE_FWK_MODULE(VertexMerger);
DEFINE_FWK_MODULE(CandidateVertexMerger);
