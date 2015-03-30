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
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

template <class VTX>
class TemplatedVertexMerger : public edm::stream::EDProducer<> {
    public:
	typedef std::vector<VTX> Product;
	TemplatedVertexMerger(const edm::ParameterSet &params);

	virtual void produce(edm::Event &event, const edm::EventSetup &es) override;

    private:
       bool trackFilter(const reco::TrackRef &track) const;
  
       void svMerger(typename Product::iterator sv, typename Product::iterator sv2);
	
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
	token_secondaryVertex = consumes<Product>(params.getParameter<edm::InputTag>("secondaryVertices"));
	produces<Product>();
}

template <class VTX>
void TemplatedVertexMerger<VTX>::produce(edm::Event &event, const edm::EventSetup &es)
{
	using namespace reco;

	edm::Handle<Product> secondaryVertices;
	event.getByToken(token_secondaryVertex, secondaryVertices);

        VertexDistance3D dist;
	std::auto_ptr<Product> recoVertices(new Product);
	for(typename Product::const_iterator sv = secondaryVertices->begin(); sv != secondaryVertices->end(); ++sv) {
          recoVertices->push_back(*sv);
        }
	
	//std::cout << "=====================================================================================================================================" << std::endl;
	//std::cout << " Original RecoVertices size is: " << recoVertices->size() <<   std::endl;


	for(typename Product::iterator sv = recoVertices->begin();  sv != recoVertices->end(); ++sv) {
	  
	  //bool shared=false;
	  VertexState s1(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error()));
	  for(typename Product::iterator sv2 = recoVertices->begin(); sv2 != recoVertices->end(); ++sv2) {
	    VertexState s2(RecoVertex::convertPos(sv2->position()),RecoVertex::convertError(sv2->error()));
	    double fr_21=vertexTools::computeSharedTracks(*sv2, *sv);
	    double fr_12=vertexTools::computeSharedTracks(*sv, *sv2);
	    // std::cout << sv2-recoVertices->begin() << " vs " << sv-recoVertices->begin() << " : " << fr << " "  <<  computeSharedTracks(*sv, *sv2) << " sig " << dist.distance(*sv,*sv2).significance() << std::endl;
	    //std::cout << "Sig all " << dist.distance(s1,s2).significance() << "Shared Tracks" << fr << std::endl;
	    
	    if (dist.distance(s1,s2).significance() < minSignificance && fr_21 > maxFraction && sv-sv2!=0 
		&& fr_21 >= fr_12){
	      //shared = true;
	      
	      //std::cout << "Sig inside " << dist.distance(s1,s2).significance() << std::endl;
	      //std::cout << "Checking " << sv-recoVertices->begin() << " and "  << sv2-recoVertices->begin() <<  std::endl;
	      
	      if(doMerging){ 
		//Merging the shared vertices
		svMerger(sv,sv2);
		
		//Erasing the merged vertex
		sv2=recoVertices->erase(sv2)-1;
	      
	      }else{ 
		//Removing the second shared vertex without merging.
		sv2=recoVertices->erase(sv2)-1;
	      
	      }
	      //std::cout << "After Merging " << sv-recoVertices->begin() << " and "  << sv2-recoVertices->begin() <<  std::endl;
	    }
	  }
	  
	  // Removing one of the shared vertices instead of merging.
	  //if(doMerging == false && shared) { 
	  //  sv=recoVertices->erase(sv)-1; 
	  //}
	//std::cout << "it = " <<  sv-recoVertices->begin() << " new size is: " << recoVertices->size() <<   std::endl;
	}
	
	//std::cout << " Final RecoVertices size is: " << recoVertices->size() <<   std::endl;
	
	event.put(recoVertices);

} //End produce


// 
template <>
void TemplatedVertexMerger<reco::Vertex>::svMerger(typename Product::iterator sv, typename Product::iterator sv2)
{
  bool sharingTracks = false;
  for(reco::Vertex::trackRef_iterator ti_sv2 = sv2->tracks_begin(); ti_sv2!= sv2->tracks_end(); ++ti_sv2){
    reco::Vertex::trackRef_iterator it = find(sv->tracks_begin(), sv->tracks_end(), *ti_sv2);
    if (it==sv->tracks_end()){
      sv->add( *ti_sv2, sv2->refittedTrack(*ti_sv2), sv2->trackWeight(*ti_sv2));
    }else sharingTracks=true;
  }
  
  if(sharingTracks){
    // create backup track containers from the main vertex
    std::vector<reco::TrackBaseRef> tracks_;
    std::vector<reco::Track> refittedTracks_;
    std::vector<float> weights_;
    
    for(reco::Vertex::trackRef_iterator it = sv->tracks_begin(); it!=sv->tracks_end(); it++) {
      tracks_.push_back( *it);
      refittedTracks_.push_back( sv->refittedTrack(*it));
      weights_.push_back( sv->trackWeight(*it) );
    }
    
    // delete tracks and add all tracks back, and check in which vertex the weight is larger
    sv->removeTracks();
    std::vector<reco::Track>::iterator it2 = refittedTracks_.begin();
    std::vector<float>::iterator it3 = weights_.begin();
    for(reco::Vertex::trackRef_iterator it = tracks_.begin(); it!=tracks_.end(); it++, it2++, it3++){
      float weight = *it3;
      float weight2= sv2->trackWeight(*it);
      reco::Track refittedTrackWithLargerWeight = *it2;
      if( weight2 >weight) {
	weight = weight2;
	refittedTrackWithLargerWeight = sv2->refittedTrack(*it);
      }
      sv->add(*it , refittedTrackWithLargerWeight , weight);
    }
  }
  
}

template <>
void TemplatedVertexMerger<reco::VertexCompositePtrCandidate>::svMerger(typename Product::iterator sv, typename Product::iterator sv2)
{
  for(std::vector<reco::CandidatePtr>::const_iterator ti_sv2 = sv2->daughterPtrVector().begin(); ti_sv2 != sv2->daughterPtrVector().end(); ++ti_sv2){
    std::vector<reco::CandidatePtr>::const_iterator it = find(sv->daughterPtrVector().begin(), sv->daughterPtrVector().end(), *ti_sv2);
    if (it==sv->daughterPtrVector().end()){ 
      sv->addDaughter(*ti_sv2);
    }
  }
}


typedef TemplatedVertexMerger<reco::Vertex> VertexMerger;
typedef TemplatedVertexMerger<reco::VertexCompositePtrCandidate> CandidateVertexMerger;

DEFINE_FWK_MODULE(VertexMerger);
DEFINE_FWK_MODULE(CandidateVertexMerger);
