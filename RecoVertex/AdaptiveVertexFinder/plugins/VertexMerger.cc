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
	//std::cout << "minSignif " << minSignificance << std::endl;
	//std::cout << "maxFraction: "<<maxFraction<<  std::endl;
	//std::cout << "Do merging: "<<doMerging<<  std::endl;

	for(typename Product::iterator sv = recoVertices->begin();  sv != recoVertices->end(); ++sv) {
	  bool mergeable=false;
	  
	  VertexState s1(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error()));
	  for(typename Product::iterator sv2 = recoVertices->begin(); sv2 != recoVertices->end(); ++sv2) {
	    VertexState s2(RecoVertex::convertPos(sv2->position()),RecoVertex::convertError(sv2->error()));
	    double fr=vertexTools::computeSharedTracks(*sv2, *sv);
	    //double fr_12=vertexTools::computeSharedTracks(*sv, *sv2);
	    // std::cout << sv2-recoVertices->begin() << " vs " << sv-recoVertices->begin() << " : " << fr << " "  <<  computeSharedTracks(*sv, *sv2) << " sig " << dist.distance(*sv,*sv2).significance() << std::endl;
	    //std::cout << "Sig all " << dist.distance(s1,s2).significance() << "Shared Tracks" << fr << std::endl;
	   
	    if (fr > maxFraction && dist.distance(s1,s2).significance() < minSignificance && sv-sv2!=0 
		&& fr >= vertexTools::computeSharedTracks(*sv, *sv2))
	      { 
	      
		mergeable= true;
		//std::cout << "In fr: "<<fr<<  std::endl;
		//std::cout << "Sig inside " << dist.distance(s1,s2).significance() << std::endl;
		//std::cout << "Checking " << sv-recoVertices->begin() << " and "  << sv2-recoVertices->begin() <<  std::endl;
		
		if(doMerging){ 
		  //Merging the shared vertices
		  svMerger(sv2,sv);
		}
	      }
	  }
	  
	  // Removing the mergeable vertex (Merged or not)
	  if(mergeable) { sv=recoVertices->erase(sv)-1; }
	  //std::cout << "it = " <<  sv-recoVertices->begin() << " new size is: " << recoVertices->size() <<   std::endl;
	}
	
	//std::cout << " Final RecoVertices size is: " << recoVertices->size() <<   std::endl;
	event.put(recoVertices);

} //End produce


// 
template <>
void TemplatedVertexMerger<reco::Vertex>::svMerger(typename Product::iterator sv2, typename Product::iterator sv)
{
  bool sharingTracks = false;
  for(reco::Vertex::trackRef_iterator ti_sv = sv->tracks_begin(); ti_sv!= sv->tracks_end(); ++ti_sv){
    reco::Vertex::trackRef_iterator it = find(sv2->tracks_begin(), sv2->tracks_end(), *ti_sv);
    if (it==sv2->tracks_end()){
      sv2->add( *ti_sv, sv->refittedTrack(*ti_sv), sv->trackWeight(*ti_sv));
    }else sharingTracks=true;
  }
  
  if(sharingTracks){
    // create backup track containers from the main vertex
    std::vector<reco::TrackBaseRef> tracks_;
    std::vector<reco::Track> refittedTracks_;
    std::vector<float> weights_;
    
    for(reco::Vertex::trackRef_iterator it = sv2->tracks_begin(); it!=sv2->tracks_end(); it++) {
      tracks_.push_back( *it);
      refittedTracks_.push_back( sv2->refittedTrack(*it));
      weights_.push_back( sv2->trackWeight(*it) );
    }
    
    // delete tracks and add all tracks back, and check in which vertex the weight is larger
    sv2->removeTracks();
    std::vector<reco::Track>::iterator it2 = refittedTracks_.begin();
    std::vector<float>::iterator it3 = weights_.begin();
    for(reco::Vertex::trackRef_iterator it = tracks_.begin(); it!=tracks_.end(); it++, it2++, it3++){
      float weight = *it3;
      float weight2= sv->trackWeight(*it);
      reco::Track refittedTrackWithLargerWeight = *it2;
      if( weight2 >weight) {
	weight = weight2;
	refittedTrackWithLargerWeight = sv->refittedTrack(*it);
      }
      sv2->add(*it , refittedTrackWithLargerWeight , weight);
    }
  }
  
}

template <>
void TemplatedVertexMerger<reco::VertexCompositePtrCandidate>::svMerger(typename Product::iterator sv2, typename Product::iterator sv)
{
  for(std::vector<reco::CandidatePtr>::const_iterator ti_sv = sv->daughterPtrVector().begin(); ti_sv != sv->daughterPtrVector().end(); ++ti_sv){
    std::vector<reco::CandidatePtr>::const_iterator it = find(sv2->daughterPtrVector().begin(), sv2->daughterPtrVector().end(), *ti_sv);
    if (it==sv2->daughterPtrVector().end()){ 
      sv2->addDaughter(*ti_sv);
    }
  }
}


typedef TemplatedVertexMerger<reco::Vertex> VertexMerger;
typedef TemplatedVertexMerger<reco::VertexCompositePtrCandidate> CandidateVertexMerger;

DEFINE_FWK_MODULE(VertexMerger);
DEFINE_FWK_MODULE(CandidateVertexMerger);
